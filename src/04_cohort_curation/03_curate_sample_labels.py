"""Curate the labels of every sample from its recorded fields by the written rules of its study platform, keeping the wording behind each value."""
import collections
import importlib
import re

C = importlib.import_module("01_catalogue_sample_metadata")

CATEGORICAL = ("sample_class", "tissue", "sex", "outcome", "septic_shock")
VARIABLES = (*CATEGORICAL, "age", "patient", "time", "phase", "definition", "scores")
# The vocabulary of a variable, where it is named otherwise than the variable.
VOCABULARY = {"age": "age_group", "phase": "sepsis_phase", "definition": "sepsis_definition"}
LABEL_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "sample_class", "infection_sepsis_not_stated", "authors_label", "class_basis",
                 "sepsis_phase", "sepsis_definition", "tissue", "sex", "age_years", "age_group", "patient_id", "time_point", "time_order", "time_unit", "outcome",
                 "outcome_definition", "septic_shock"]
# The phases that a sample of each kind may have: a sepsis sample is drawn before, during or at an unrecorded time relative to the onset; a follow-up or
# recovery sample of a sepsis patient is of another condition.
SEPSIS_PHASES = ("before_onset", "during_sepsis", "onset_not_recorded")
OTHER_PHASES = ("recovery_or_follow_up", "not_applicable")
SOURCE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "variable", "value", "basis", "rule", "fields", "recorded", "matched",
                  "study_statement", "note"]
SCORE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "score", "value", "field", "recorded"]
REVIEW_COLUMNS = ["study_accession", "platform", "variable", "rule", "fields", "matched", "value", "samples", "lowest", "highest", "examples", "basis",
                  "study_statement", "note"]
NOT_RECORDED = "not_recorded"
# Days in one unit of a recorded age: exact in binary, so that an age recorded at a limit is grouped at that limit.
AGE_DAYS = {"years": 365.25, "months": 365.25 / 12, "weeks": 7.0, "days": 1.0}


def plain(text):
    """Write a recorded wording in lower case with single spaces, the form in which a rule names it."""
    return re.sub(r"\s+", " ", str(text)).strip().lower()


def rule_list(entry):
    return entry if isinstance(entry, list) else [entry]


def field_list(rule):
    names = rule.get("from", [])
    return [names] if isinstance(names, str) else list(names)


def wording(fields, names):
    """Return what one sample records in the named fields, joined in their order, or None where it records none of them."""
    values = [fields.get(name) for name in names]
    return None if all(value is None for value in values) else " | ".join(value or "" for value in values)


def decide(rule, text, shared):
    """Apply one rule to a recorded wording and return the key that matched and the value it gives, or None where the rule leaves the sample to the next rule."""
    if "values" in rule:
        mapping = shared[rule["values"]] if isinstance(rule["values"], str) else rule["values"]
        key = plain(text)
        if key not in mapping:
            C.require("other" in rule, f"A recorded value is not named in its rule: {text!r} in {field_list(rule)}")
            return "(any other value)", rule["other"]
        return key, mapping[key]
    if "patterns" in rule:
        for pattern, value in rule["patterns"]:
            if re.search(pattern, text, re.IGNORECASE):
                return pattern, value
        return None
    pattern = rule["number"] if "number" in rule else rule["text"]
    match = re.search(pattern, text, re.IGNORECASE)
    if not match:
        return None
    return pattern, float(match.group(1)) if "number" in rule else "_".join(group for group in match.groups() if group)


def resolve(entry, fields, shared):
    """Return the decision of the first rule that gives one sample a value: rule number, wording, matched key, value and the rule itself."""
    for number, rule in enumerate(rule_list(entry), 1):
        names = field_list(rule)
        if not names:
            return number, "", "(study description)", rule["value"], rule
        text = wording(fields, names)
        if text is None:
            continue
        decision = decide(rule, text, shared)
        if decision is not None and decision[1] != "next":
            return number, text, decision[0], decision[1], rule
    return 0, "", "", None, {}


def basis_of(rule, number):
    if not number:
        return NOT_RECORDED
    if not field_list(rule):
        return "study_description"
    return "sample_record_and_study_description" if "study_states" in rule else "sample_record"


def age_group(days, limits):
    """Group an age from its unrounded number of days."""
    if days <= limits["neonate_up_to_days"]:
        return "neonate"
    return "adult" if days >= limits["adult_from_years"] * AGE_DAYS["years"] else "child"


def check_rules(key, rules, present, study_text, policy):
    """Stop where a rule set names an unknown variable, a field its study platform does not record, a value outside the vocabulary or a sentence its study does not hold."""
    vocabularies = policy["vocabularies"]
    C.require(set(rules) <= set(VARIABLES) - {"definition"}, f"Unknown variable in the rules of {key}: {sorted(set(rules) - set(VARIABLES))}")
    entries = [(variable, rule) for variable, entry in rules.items() if variable != "scores" for rule in rule_list(entry)]
    entries += [("scores", rule) for entry in rules.get("scores", {}).values() for rule in rule_list(entry)]
    for variable, rule in entries:
        names = field_list(rule)
        C.require(set(names) <= present, f"A rule of {key} names a field that its samples do not record: {sorted(set(names) - present)}")
        C.require(sum(name in rule for name in ("values", "patterns", "number", "text")) == (1 if names else 0) and (names or {"study_states", "value"} <= set(rule)),
                  f"A rule of {key} for {variable} must read fields in one way, or state a value from the study description")
        if "study_states" in rule:
            C.require(plain(rule["study_states"]) in study_text, f"The study record of {key} does not hold the sentence: {rule['study_states']!r}")
        if variable in vocabularies or variable == "phase" or (variable == "age" and "number" not in rule):
            allowed = set(vocabularies[VOCABULARY.get(variable, variable)]) | {"next", None}
            mapping = policy["shared"][rule["values"]] if isinstance(rule.get("values"), str) else rule.get("values", {})
            named = [*mapping.values(), *(value for _, value in rule.get("patterns", [])), *([rule["other"]] if "other" in rule else []),
                     *([rule["value"]] if "value" in rule else [])]
            C.require(set(named) <= allowed, f"A rule of {key} for {variable} gives a value outside the vocabulary: {sorted(set(map(str, named)) - allowed)}")
        if variable == "age" and "number" in rule:
            C.require(rule.get("unit") in AGE_DAYS, f"An age rule of {key} needs a unit")
        if variable == "time":
            C.require(rule.get("unit") in vocabularies["time_unit"], f"A time rule of {key} needs a unit")


def main():
    import numpy as np
    import pandas as pd

    print(f"[03] Curating sample labels v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy = C.CONFIG["curation"]
    shared, limits = policy["shared"], policy["age_groups"]
    metadata = pd.read_parquet(C.DATA / "sample_metadata.parquet")
    index = C.sample_records(["source", "study_accession", "sample_accession", "platform"])
    studies = pd.read_parquet(C.ACQUIRED / "studies.parquet", columns=["accession", "title", "summary", "design"]).set_index("accession")
    tables = sorted(set(zip(index["study_accession"], index["platform"])))
    keys = {f"{study}|{platform}" for study, platform in tables}
    C.require(set(C.REGISTER["rules"]) == keys, f"The rules of {C.CONFIG['register']} must name every study platform with expression values and no "
                                                f"other: {sorted(set(C.REGISTER['rules']) ^ keys)}")
    definitions = C.REGISTER["sepsis_definitions"]
    C.require(set(definitions) <= {study for study, _ in tables} and all(entry["value"] in policy["vocabularies"]["sepsis_definition"] for entry in definitions.values()),
              "A recorded sepsis definition names a study that is not included, or a value outside its vocabulary")

    # What every sample records in every field: the values of one field in their recorded order, each once.
    recorded = collections.defaultdict(lambda: collections.defaultdict(dict))
    # Of the fields of the record itself only the free-text ones are read, so that none is confused with a characteristic of the same name.
    metadata = metadata[(metadata["origin"] != "record") | metadata["field"].isin(policy["record_fields"])]
    for study, platform, sample, field, value in metadata[["study_accession", "platform", "sample_accession", "field", "value"]].itertuples(index=False, name=None):
        recorded[(study, platform)][sample].setdefault(field, {})[str(value).strip()] = None
    labels, sources, scores, used = [], [], [], collections.Counter()
    for study, platform in tables:
        key = f"{study}|{platform}"
        rules = C.REGISTER["rules"][key]
        part = index[(index["study_accession"] == study) & (index["platform"] == platform)]
        samples = {sample: {field: " + ".join(values) for field, values in fields.items()} for sample, fields in recorded[(study, platform)].items()}
        C.require(set(part["sample_accession"]) == set(samples), f"The catalogue does not hold every sample of {key}")
        row = studies.loc[study]
        study_text = plain(f"{row['title']}. {row['summary']} {row['design']}")
        check_rules(key, rules, {field for fields in samples.values() for field in fields}, study_text, policy)
        definition = definitions.get(study)
        C.require(definition is None or plain(definition["study_states"]) in study_text, f"The study record of {study} does not hold the sentence of its sepsis definition")
        for item in part.itertuples(index=False):
            fields = samples[item.sample_accession]
            label = dict(source=item.source, study_accession=study, sample_accession=item.sample_accession, platform=platform, sample_class=NOT_RECORDED,
                         infection_sepsis_not_stated=False, authors_label="", class_basis=NOT_RECORDED, sepsis_phase="", sepsis_definition="",
                         tissue=NOT_RECORDED, sex=NOT_RECORDED, age_years=np.nan,
                         age_group=NOT_RECORDED, patient_id="", time_point="", time_order=np.nan, time_unit="", outcome=NOT_RECORDED, outcome_definition="",
                         septic_shock=NOT_RECORDED)
            for variable in VARIABLES[:-1]:
                if variable == "definition":
                    # The definition of sepsis is a property of the study, recorded only where its record states it.
                    number, text, matched, rule = (1, "", "(study description)", definition) if definition else (0, "", "", {})
                    value = definition["value"] if definition else "not_stated_in_record"
                else:
                    number, text, matched, value, rule = resolve(rules.get(variable, []), fields, shared)
                    used[(key, variable, number, matched)] += 1
                basis = "study_description" if variable == "definition" and number else basis_of(rule, number)
                if variable == "phase" and value is None:
                    # A sample whose rule set records no phase is during sepsis if it is a sepsis sample, and outside it otherwise.
                    value, basis = ("during_sepsis" if label["sample_class"] == "sepsis" else "not_applicable"), "from_sample_class"
                value = NOT_RECORDED if value is None else value
                if variable == "sample_class":
                    label.update(sample_class="other_condition" if value == "infection_sepsis_not_stated" else value, authors_label=text,
                                 infection_sepsis_not_stated=value == "infection_sepsis_not_stated", class_basis=basis_of(rule, number))
                elif variable == "age" and isinstance(value, float):
                    days = value * AGE_DAYS[rule["unit"]]
                    label.update(age_years=round(days / AGE_DAYS["years"], 4), age_group=age_group(days, limits))
                elif variable == "age":
                    label.update(age_group=value)
                elif variable == "patient":
                    label.update(patient_id="" if value == NOT_RECORDED else value)
                elif variable == "time" and value != NOT_RECORDED:
                    label.update(time_point=text, time_order=float(value), time_unit=rule["unit"])
                elif variable == "outcome":
                    label.update(outcome=value, outcome_definition=rule.get("definition", "") if value != NOT_RECORDED else "")
                elif variable in ("tissue", "sex", "septic_shock"):
                    label[variable] = value
                elif variable == "phase":
                    C.require(value in (SEPSIS_PHASES if label["sample_class"] == "sepsis" else OTHER_PHASES)
                              and (value != "recovery_or_follow_up" or label["sample_class"] == "other_condition"),
                              f"{key}: sample {item.sample_accession} of class {label['sample_class']} cannot have the phase {value}")
                    label["sepsis_phase"] = value
                elif variable == "definition":
                    label["sepsis_definition"] = value
                sources.append(dict(source=item.source, study_accession=study, sample_accession=item.sample_accession, platform=platform, variable=variable,
                                    value=str(value), basis=basis, rule=number, fields=" | ".join(field_list(rule)), recorded=text, matched=matched,
                                    study_statement=rule.get("study_states", ""), note=rule.get("note", "")))
            for score, entry in rules.get("scores", {}).items():
                number, text, matched, value, rule = resolve(entry, fields, shared)
                used[(key, f"scores.{score}", number, matched)] += 1
                if value is not None:
                    scores.append(dict(source=item.source, study_accession=study, sample_accession=item.sample_accession, platform=platform, score=score,
                                       value=float(value), field=" | ".join(field_list(rule)), recorded=text))
            labels.append(label)
        # A rule that no sample of its study platform meets is a mistake in the rule, not a finding.
        for variable, entry in [(name, rules[name]) for name in rules if name != "scores"] + [(f"scores.{name}", rule) for name, rule in rules.get("scores", {}).items()]:
            for number, rule in enumerate(rule_list(entry), 1):
                named = rule["values"].items() if isinstance(rule.get("values"), dict) else rule.get("patterns", [])
                named = [name for name, value in named if value != "next"]
                idle = [name for name in named if not used[(key, variable, number, name)]]
                met = any(count for (table, name, rule_number, _), count in used.items() if (table, name, rule_number) == (key, variable, number))
                C.require(not idle and (met or not named and "values" in rule), f"A rule of {key} for {variable} is met by no sample: rule {number} {idle}")

    labels = pd.DataFrame(labels, columns=LABEL_COLUMNS)
    sources = pd.DataFrame(sources, columns=SOURCE_COLUMNS)
    scores = pd.DataFrame(scores, columns=SCORE_COLUMNS)
    C.require(len(labels) == len(index) and set(labels["sample_class"]) <= set(policy["vocabularies"]["sample_class"]), "Every sample must have one curated class")
    numeric = pd.to_numeric(sources["value"], errors="coerce").where(sources["variable"].isin(["age", "time"]))
    review = (sources.assign(number=numeric, value=sources["value"].where(numeric.isna(), "(number)").where(sources["variable"] != "patient", "(identifier)"))
              .groupby(["study_accession", "platform", "variable", "rule", "fields", "matched", "value", "basis", "study_statement", "note"], sort=False, dropna=False)
              .agg(samples=("sample_accession", "size"), lowest=("number", "min"), highest=("number", "max"),
                   examples=("recorded", lambda values: "; ".join(str(value)[:60] for value in list(dict.fromkeys(values))[:3]))).reset_index())
    review = review.sort_values(["study_accession", "platform", "variable", "rule", "matched"], kind="stable")[REVIEW_COLUMNS]
    C.save_table(labels, C.DATA / "sample_labels.parquet")
    C.save_table(sources, C.DATA / "label_sources.parquet")
    C.save_table(scores, C.DATA / "sample_scores.parquet")
    C.save_review(review, C.OUT / "label_rules_review.csv")

    def counts(column):
        return {str(name): int(count) for name, count in labels[column].value_counts().sort_index().items()}

    by_table = labels.groupby(["study_accession", "platform"])
    C.finish(3, start, [C.DATA / "sample_labels.parquet", C.DATA / "label_sources.parquet", C.DATA / "sample_scores.parquet", C.OUT / "label_rules_review.csv"],
             method=policy["method"], samples=len(labels), study_platforms=len(tables), rule_sets=len(C.REGISTER["rules"]),
             sample_class=counts("sample_class"), class_basis=counts("class_basis"), infection_sepsis_not_stated=int(labels["infection_sepsis_not_stated"].sum()),
             sepsis_phase=counts("sepsis_phase"), sepsis_definition=counts("sepsis_definition"),
             sepsis_by_phase={str(name): int(count) for name, count in labels.loc[labels["sample_class"] == "sepsis", "sepsis_phase"].value_counts().sort_index().items()},
             tissue=counts("tissue"), sex=counts("sex"), age_group=counts("age_group"), outcome=counts("outcome"), septic_shock=counts("septic_shock"),
             samples_with=dict(age_in_years=int(labels["age_years"].notna().sum()), patient_identifier=int((labels["patient_id"] != "").sum()),
                               time_point=int(labels["time_order"].notna().sum()), score=int(scores["sample_accession"].nunique()) if len(scores) else 0),
             study_platforms_with=dict(sepsis=int((by_table["sample_class"].agg(lambda values: (values == "sepsis").any())).sum()),
                                       healthy_control=int((by_table["sample_class"].agg(lambda values: (values == "healthy_control").any())).sum()),
                                       patient_identifier=int((by_table["patient_id"].agg(lambda values: (values != "").any())).sum()),
                                       time_point=int((by_table["time_order"].agg(lambda values: values.notna().any())).sum())),
             scores={str(name): int(count) for name, count in scores["score"].value_counts().sort_index().items()},
             vocabularies=policy["vocabularies"], age_groups=limits, protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: {len(labels):,} samples curated in {len(tables)} study platforms; class: "
          + ", ".join(f"{name} {count:,}" for name, count in counts("sample_class").items())
          + f"; review sheet: {(C.OUT / 'label_rules_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
