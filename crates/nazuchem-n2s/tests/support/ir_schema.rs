// SPDX-License-Identifier: MIT
// Origin: original
// Source: Nazuchem | crates/nazuchem-n2s/tests/support/ir_schema.rs | Phase B step 8a

//! Offline validation of the published IR JSON Schema, independent of serde.
//!
//! This deliberately implements only the schema vocabulary exercised by the
//! published contract. Preflight rejects any unknown keyword, including one in
//! an unused definition, so a future schema cannot silently weaken this gate.
//! Annotation keywords (including `format`, as in the default 2020-12 dialect)
//! impose no additional assertions. Only local JSON Pointer references are used.

use serde_json::Value;

pub struct Validator {
    schema: Value,
}

impl Validator {
    pub fn new(schema: Value) -> Result<Self, String> {
        let validator = Self { schema };
        validator.preflight(&validator.schema, "schema")?;
        Ok(validator)
    }

    pub fn validate(&self, instance: &Value) -> Result<(), String> {
        self.check(&self.schema, instance, "$", 0)
    }

    fn reference(&self, reference: &str) -> Result<&Value, String> {
        if reference == "#" {
            return Ok(&self.schema);
        }
        let pointer = reference
            .strip_prefix('#')
            .filter(|pointer| pointer.starts_with('/'))
            .ok_or_else(|| format!("unsupported nonlocal schema reference {reference}"))?;
        self.schema
            .pointer(pointer)
            .ok_or_else(|| format!("unresolved schema reference {reference}"))
    }

    fn preflight(&self, schema: &Value, path: &str) -> Result<(), String> {
        if schema.is_boolean() {
            return Ok(());
        }
        let object = schema
            .as_object()
            .ok_or_else(|| format!("{path}: schema must be an object or boolean"))?;
        for (keyword, value) in object {
            let location = format!("{path}.{keyword}");
            let valid = match keyword.as_str() {
                "$schema" | "$id" | "title" | "description" | "format" => value.is_string(),
                "default" | "const" => true,
                "examples" | "enum" => value.is_array(),
                "deprecated" | "readOnly" | "writeOnly" | "uniqueItems" => value.is_boolean(),
                "$ref" => {
                    let reference = value
                        .as_str()
                        .ok_or_else(|| format!("{location}: expected a string"))?;
                    self.reference(reference)?;
                    true
                }
                "$defs" | "properties" => {
                    let definitions = value
                        .as_object()
                        .ok_or_else(|| format!("{location}: expected an object"))?;
                    for (name, child) in definitions {
                        self.preflight(child, &format!("{location}.{name}"))?;
                    }
                    true
                }
                "items" | "additionalProperties" | "not" => {
                    self.preflight(value, &location)?;
                    true
                }
                "allOf" | "anyOf" | "oneOf" | "prefixItems" => {
                    let choices = value
                        .as_array()
                        .ok_or_else(|| format!("{location}: expected an array"))?;
                    for (index, child) in choices.iter().enumerate() {
                        self.preflight(child, &format!("{location}[{index}]"))?;
                    }
                    true
                }
                "type" => {
                    let supported = |kind: &Value| {
                        matches!(
                            kind.as_str(),
                            Some(
                                "null"
                                    | "boolean"
                                    | "object"
                                    | "array"
                                    | "number"
                                    | "integer"
                                    | "string"
                            )
                        )
                    };
                    supported(value)
                        || value
                            .as_array()
                            .is_some_and(|types| !types.is_empty() && types.iter().all(supported))
                }
                "required" => value
                    .as_array()
                    .is_some_and(|names| names.iter().all(Value::is_string)),
                "minItems" | "maxItems" | "minLength" | "maxLength" | "minProperties"
                | "maxProperties" => value.as_u64().is_some(),
                "minimum" | "maximum" | "exclusiveMinimum" | "exclusiveMaximum" => {
                    value.is_number()
                }
                _ => return Err(format!("{location}: unsupported JSON Schema keyword")),
            };
            if !valid {
                return Err(format!("{location}: invalid keyword value"));
            }
        }
        Ok(())
    }

    fn check(
        &self,
        schema: &Value,
        instance: &Value,
        path: &str,
        depth: usize,
    ) -> Result<(), String> {
        if depth > 256 {
            return Err(format!("{path}: schema recursion limit"));
        }
        match schema.as_bool() {
            Some(true) => return Ok(()),
            Some(false) => return Err(format!("{path}: forbidden by false schema")),
            None => {}
        }
        let object = schema.as_object().expect("preflight checked schema shape");
        if let Some(reference) = object.get("$ref") {
            self.check(
                self.reference(reference.as_str().unwrap())?,
                instance,
                path,
                depth + 1,
            )?;
        }
        let fail = |keyword: &str| Err(format!("{path}: violates {keyword}"));
        if let Some(kinds) = object.get("type") {
            let matches = |kind: &Value| match kind.as_str().unwrap() {
                "null" => instance.is_null(),
                "boolean" => instance.is_boolean(),
                "object" => instance.is_object(),
                "array" => instance.is_array(),
                "number" => instance.is_number(),
                "integer" => instance
                    .as_f64()
                    .is_some_and(|number| number.fract() == 0.0),
                "string" => instance.is_string(),
                _ => unreachable!("preflight checked types"),
            };
            if !kinds
                .as_array()
                .map_or_else(|| matches(kinds), |kinds| kinds.iter().any(matches))
            {
                return fail("type");
            }
        }
        if object.get("const").is_some_and(|value| instance != value) {
            return fail("const");
        }
        if object
            .get("enum")
            .is_some_and(|values| !values.as_array().unwrap().contains(instance))
        {
            return fail("enum");
        }
        for keyword in ["allOf", "anyOf", "oneOf"] {
            if let Some(choices) = object.get(keyword) {
                let choices = choices.as_array().unwrap();
                let valid = choices
                    .iter()
                    .filter(|child| self.check(child, instance, path, depth + 1).is_ok())
                    .count();
                if !match keyword {
                    "allOf" => valid == choices.len(),
                    "anyOf" => valid > 0,
                    _ => valid == 1,
                } {
                    return fail(keyword);
                }
            }
        }
        if object
            .get("not")
            .is_some_and(|child| self.check(child, instance, path, depth + 1).is_ok())
        {
            return fail("not");
        }
        if let Some(properties) = instance.as_object() {
            if let Some(required) = object.get("required") {
                for name in required.as_array().unwrap() {
                    let name = name.as_str().unwrap();
                    if !properties.contains_key(name) {
                        return Err(format!("{path}.{name}: missing required property"));
                    }
                }
            }
            for (name, value) in properties {
                if let Some(child) = object.get("properties").and_then(|p| p.get(name)) {
                    self.check(child, value, &format!("{path}.{name}"), depth + 1)?;
                } else if let Some(child) = object.get("additionalProperties") {
                    self.check(child, value, &format!("{path}.{name}"), depth + 1)?;
                }
            }
            for (keyword, minimum) in [("minProperties", true), ("maxProperties", false)] {
                if object.get(keyword).is_some_and(|limit| {
                    if minimum {
                        properties.len() < limit.as_u64().unwrap() as usize
                    } else {
                        properties.len() > limit.as_u64().unwrap() as usize
                    }
                }) {
                    return fail(keyword);
                }
            }
        }
        if let Some(items) = instance.as_array() {
            let prefix = object.get("prefixItems").and_then(Value::as_array);
            for (index, item) in items.iter().enumerate() {
                let child = prefix
                    .and_then(|prefix| prefix.get(index))
                    .or_else(|| object.get("items"));
                if let Some(child) = child {
                    self.check(child, item, &format!("{path}[{index}]"), depth + 1)?;
                }
            }
            for (keyword, minimum) in [("minItems", true), ("maxItems", false)] {
                if object.get(keyword).is_some_and(|limit| {
                    if minimum {
                        items.len() < limit.as_u64().unwrap() as usize
                    } else {
                        items.len() > limit.as_u64().unwrap() as usize
                    }
                }) {
                    return fail(keyword);
                }
            }
            if object.get("uniqueItems") == Some(&Value::Bool(true))
                && items
                    .iter()
                    .enumerate()
                    .any(|(index, value)| items[..index].contains(value))
            {
                return fail("uniqueItems");
            }
        }
        if let Some(string) = instance.as_str() {
            let length = string.chars().count();
            for (keyword, minimum) in [("minLength", true), ("maxLength", false)] {
                if object.get(keyword).is_some_and(|limit| {
                    if minimum {
                        length < limit.as_u64().unwrap() as usize
                    } else {
                        length > limit.as_u64().unwrap() as usize
                    }
                }) {
                    return fail(keyword);
                }
            }
        }
        if let Some(number) = instance.as_f64() {
            for keyword in ["minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum"] {
                if let Some(limit) = object.get(keyword).and_then(Value::as_f64)
                    && match keyword {
                        "minimum" => number < limit,
                        "maximum" => number > limit,
                        "exclusiveMinimum" => number <= limit,
                        _ => number >= limit,
                    }
                {
                    return fail(keyword);
                }
            }
        }
        Ok(())
    }
}
