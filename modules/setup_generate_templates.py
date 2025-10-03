#!/usr/bin/env python3
"""
Creates:
  1) hmm_table_template.txt
  2) hmm_table_template_2.txt  (collapsed; merges KO/HMM/threshold/abbrs)
  3) kegg_module_step_db.txt   (module steps as boolean KO expressions)

Inputs (adjust paths below if needed):
  - databases/curated.txt
      columns: Category, EC_or_KO, gene_name(s ; separated), hmm_file, hmm_threshold
  - databases/kofam_database/ko_list
      tab-delimited: KO<TAB>threshold (numeric or pattern), one per line
  - databases/kofam_database/profiles/prokaryote.hal
      list of KO *.hmm files shipped with KOfam (one per line)
  - templates/ko00002.json

Notes
-----
- Collapsing uses (Function, Category) as the grouping key, preserving first-seen order.
- When collapsing, KO/HMM/Threshold/Gene_abbreviation are **merged (unique)** and comma-joined.
- KO thresholds prefer the KO-level threshold from ko_list when present; otherwise fallback to DEFAULT_THRESHOLD.
"""

import os
import re
import json
import time
import requests
import boolean  # pip install boolean.py
from pathlib import Path
from collections import defaultdict, OrderedDict
from tqdm import tqdm

# -----------------------
# Paths / constants
# -----------------------

INPUT_FILE = "databases/curated.txt"
KO_LIST_FILE = "databases/kofam_database/ko_list"
OUTPUT_FILE = "templates/hmm_table_template.txt"
OUTPUT_COLLAPSED = "templates/hmm_table_template_2.txt"
DEFAULT_THRESHOLD = "000.00|full"

# For module step DB
HAL_PATH = "databases/kofam_database/profiles/prokaryote.hal"
KO_JSON = "templates/ko00002.json"          # pathway tree (used to list modules)
OUTPUT_MODULE_STEPS = "templates/kegg_module_step_db.txt"

# -----------------------
# Globals / caches
# -----------------------
reaction_cache = {}
ko_info_cache = {}
ko_thresholds = {}

# -----------------------
# Utilities: thresholds
# -----------------------
def load_ko_thresholds(path: str) -> dict:
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                d[parts[0]] = parts[1] + "|full"
    return d

ko_thresholds = load_ko_thresholds(KO_LIST_FILE)

# -----------------------
# KEGG helpers
# -----------------------
def resolve_ec(ec: str) -> str:
    """Follow 'Now EC'/'Transferred to' mapping; fallback to original EC."""
    try:
        res = requests.get(f"https://rest.kegg.jp/get/ec:{ec}", timeout=30)
        text = res.text
        match = re.search(r"(?:Now EC|Transferred to) ([\d\.]+)", text)
        if match:
            return match.group(1)
    except Exception:
        pass
    return ec

def get_kos_for_ec(ec: str):
    try:
        res = requests.get(f"https://rest.kegg.jp/link/ko/ec:{ec}", timeout=30)
        return [
            line.split("\t")[1].replace("ko:", "")
            for line in res.text.strip().split("\n")
            if "ko:" in line
        ]
    except Exception:
        return []

def get_reaction_definition(rxn_id: str):
    if rxn_id in reaction_cache:
        return reaction_cache[rxn_id]
    try:
        res = requests.get(f"https://rest.kegg.jp/get/{rxn_id}", timeout=30)
        for line in res.text.splitlines():
            if line.startswith("DEFINITION"):
                eq = line.replace("DEFINITION", "").strip()
                reaction_cache[rxn_id] = eq.replace("<=>", "=").replace("=>", "=") + f" [RN:{rxn_id}]"
                return reaction_cache[rxn_id]
    except Exception:
        pass
    return None

def get_ko_info(ko: str):
    """Return dict with Gene_name, Gene_abbreviation (list), Reaction, Substrate, Product."""
    if ko in ko_info_cache:
        return ko_info_cache[ko]
    try:
        res = requests.get(f"https://rest.kegg.jp/get/{ko}", timeout=30)
        text = res.text
        name = "Unknown"
        symbols = []
        reactions, substrates, products = [], [], []

        for line in text.splitlines():
            if line.startswith("NAME"):
                name = re.sub(r"\s*\[EC:[^\]]+\]", "", line.replace("NAME", "").strip())
            elif line.startswith("SYMBOL"):
                # avoid accidental EC-like strings
                symbols = [
                    s.strip() for s in line.replace("SYMBOL", "").split(",")
                    if not re.match(r"^E\d+\.\d+\.\d+\.\d+$", s.strip())
                ]
            elif "REACTION" in line or "EQUATION" in line:
                rxn = re.search(r"(R\d{5})", line)
                if rxn:
                    eq = get_reaction_definition(rxn.group(1))
                    if eq:
                        reactions.append(eq)
                        parts = eq.split(" [RN")[0].split("=")
                        if len(parts) == 2:
                            substrates += [s.strip() for s in parts[0].split("+")]
                            products  += [p.strip() for p in parts[1].split("+")]

        ko_info_cache[ko] = {
            "Gene_name": name,
            "Gene_abbreviation": symbols or ["Unknown"],
            "Reaction": " ; ".join(reactions) or "N/A",
            "Substrate": "; ".join(substrates) or "N/A",
            "Product": "; ".join(products) or "N/A",
        }
        time.sleep(0.2)  # helps with KEGG API
        return ko_info_cache[ko]
    except Exception:
        return {
            "Gene_name": "Unknown",
            "Gene_abbreviation": ["Unknown"],
            "Reaction": "N/A",
            "Substrate": "N/A",
            "Product": "N/A",
        }

def get_first_pathway_for_ec(ec: str) -> str:
    try:
        url = f"https://rest.kegg.jp/get/ec:{ec}"
        res = requests.get(url, timeout=30)
        if res.status_code != 200:
            return "Unknown KEGG Pathway"
        for line in res.text.splitlines():
            if line.startswith("PATHWAY"):
                # "PATHWAY mapXXXX  Pathway name"
                return line.split(None, 2)[2].strip()
    except Exception:
        pass
    return "Unknown KEGG Pathway"

# -----------------------
# Build HMM table template
# -----------------------
def build_hmm_table_template():
    if os.path.exists(OUTPUT_FILE):
        print(f"⚠️ {OUTPUT_FILE} already exists. Skipping regeneration.")
        return

    print(f"Building HMM table template: {OUTPUT_FILE}")

    SHOW_INNER_BARS = False  # set True for a nested KO bar (can be noisy)

    with open(INPUT_FILE, encoding="latin1") as f, open(OUTPUT_FILE, "w", encoding="utf-8") as out:
        out.write("#Entry\tCategory\tFunction\tEC\tGene_abbreviation\tGene_name\tKO\tReaction\tSubstrate\tProduct\tHmm_file\tHmm_detecting_threshold\n")

        # Pre-read + filter for total
        print("Counting input lines...")
        raw_lines = f.readlines()
        lines = [ln for ln in raw_lines if ln.strip() and not ln.lower().startswith("category")]
        print(f"Found {len(lines)} entries in {INPUT_FILE}")

        entry_counter = 1
        for line in tqdm(lines, desc="Building HMM table", unit="entry"):
            parts = line.strip().split("\t")
            while len(parts) < 5:
                parts.append("N/A")

            category, ec_raw, gene_name_str, hmm_file, hmm_thresh = parts[:5]
            gene_names = [g.strip().lower() for g in gene_name_str.split(";")]

            matched_kos, matched_abbrs, matched_thresholds, hmm_files = [], [], [], []
            reaction = substrate = product = "N/A"
            function = "Unknown KEGG Pathway"

            if hmm_file != "N/A":
                hmm_files.append(hmm_file)
                matched_thresholds.append(hmm_thresh if hmm_thresh and hmm_thresh != "N/A" else DEFAULT_THRESHOLD)

            kos = []
            ec_label_for_bar = None
            if ec_raw != "N/A":
                if ec_raw.startswith("K"):
                    kos = [ec_raw]
                    function = get_ko_info(ec_raw)["Gene_name"]
                else:
                    ec_resolved = resolve_ec(ec_raw)
                    kos = get_kos_for_ec(ec_resolved)
                    function = get_first_pathway_for_ec(ec_resolved)
                    ec_label_for_bar = ec_resolved

            ko_iter = kos
            if SHOW_INNER_BARS and kos:
                ko_iter = tqdm(kos, desc=f"KOs for {ec_label_for_bar or 'KO'}", unit="KO", leave=False)

            for ko in ko_iter:
                info = get_ko_info(ko)
                ko_gene_names = [k.strip().lower() for k in info["Gene_name"].split(",")]
                if any(g in ko_gene_names for g in gene_names) or not gene_names or gene_names == ["n/a"]:
                    matched_kos.append(ko)
                    matched_abbrs += info["Gene_abbreviation"]
                    ko_hmm = f"{ko}.hmm"
                    if ko_hmm not in hmm_files:
                        hmm_files.append(ko_hmm)
                        matched_thresholds.append(ko_thresholds.get(ko, DEFAULT_THRESHOLD))
                    if reaction == "N/A":
                        reaction = info["Reaction"]
                        substrate = info["Substrate"]
                        product = info["Product"]

            entry_id = f"{entry_counter:03d}"
            out.write("\t".join([
                entry_id,
                category,
                function,
                ec_raw,
                ", ".join(sorted(set([a for a in matched_abbrs if a and a != 'Unknown']))) if matched_abbrs else "Unknown",
                gene_name_str,
                ", ".join(sorted(set(matched_kos))) if matched_kos else "N/A",
                reaction,
                substrate,
                product,
                ", ".join(hmm_files) if hmm_files else "N/A",
                ", ".join(matched_thresholds) if matched_thresholds else DEFAULT_THRESHOLD
            ]) + "\n")

            entry_counter += 1

    print(f"✅ Wrote {OUTPUT_FILE}")

# -----------------------
# Collapse table (merge KO/HMM/threshold/abbrs)
# -----------------------
def collapse_hmm_table():
    category_function_order = OrderedDict()
    merged = defaultdict(lambda: {
        "entries": [], "abbrs": set(), "kos": set(), "hmms": set(),
        "thresholds": set(), "category": None, "function": None,
    })

    # Count lines for total (minus header)
    with open(OUTPUT_FILE, encoding="utf-8") as _counter:
        total_rows = max(sum(1 for _ in _counter) - 1, 0)

    with open(OUTPUT_FILE, encoding="utf-8") as infile:
        header = next(infile, None)
        row_iter = tqdm(infile, total=total_rows, desc="Collapsing rows", unit="row")
        for line in row_iter:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            entry_id = parts[0]
            category = parts[1]
            function = parts[2]
            abbrs = parts[4]
            ko_str = parts[6]
            hmm_str = parts[10]
            thr_str = parts[11]

            key = (function, category)
            if key not in category_function_order:
                category_function_order[key] = len(category_function_order)

            merged[key]["entries"].append(entry_id)
            merged[key]["category"] = category
            merged[key]["function"] = function

            if abbrs and abbrs != "Unknown":
                for a in abbrs.split(","):
                    merged[key]["abbrs"].add(a.strip())

            if ko_str and ko_str != "N/A":
                for k in ko_str.split(","):
                    k = k.strip()
                    if k:
                        merged[key]["kos"].add(k)

            if hmm_str and hmm_str != "N/A":
                for h in hmm_str.split(","):
                    h = h.strip()
                    if h:
                        merged[key]["hmms"].add(h)

            if thr_str and thr_str != "N/A":
                for t in thr_str.split(","):
                    t = t.strip()
                    if t:
                        merged[key]["thresholds"].add(t)

    with open(OUTPUT_COLLAPSED, "w", encoding="utf-8") as out2:
        out2.write("#Row\tEntry\tCategory\tFunction\tGene_abbreviation\tKO\tHmm_file\tHmm_detecting_threshold\n")
        for i, key in enumerate(category_function_order.keys(), 1):
            data = merged[key]
            out2.write("\t".join([
                f"{i:03d}",
                "||".join(data["entries"]),
                data["category"],
                data["function"],
                ", ".join(sorted({a for a in data["abbrs"] if a})),
                ", ".join(sorted({k for k in data["kos"] if k})),
                ", ".join(sorted({h for h in data["hmms"] if h})),
                ", ".join(sorted({t for t in data["thresholds"] if t})) if data["thresholds"] else DEFAULT_THRESHOLD
            ]) + "\n")

    print(f"✅ Wrote {OUTPUT_COLLAPSED}")

# -----------------------
# KEGG Module step DB
# -----------------------
algebra = boolean.BooleanAlgebra()

def load_kos_from_hal(path):
    with open(path) as f:
        return {line.strip().replace('.hmm', '') for line in f if line.startswith('K')}

def extract_modules(json_path):
    with open(json_path) as f:
        data = json.load(f)

    def collect_modules(node):
        modules = []
        if isinstance(node, dict):
            name = node.get("name", "")
            match = re.match(r"^(M\d{5})\s+(.*?)(?:\s+\[.*?\])?$", name)
            if match:
                modules.append((match.group(1), match.group(2)))
            for child in node.get("children", []):
                modules.extend(collect_modules(child))
        return modules

    return collect_modules(data)

def get_kegg_module_definition(module_id):
    url = f"http://rest.kegg.jp/get/{module_id}"
    response = requests.get(url, timeout=30)
    if response.status_code == 200:
        for line in response.text.splitlines():
            if line.startswith("DEFINITION"):
                return line.replace("DEFINITION", "").strip()
    return None

def normalize_expression(expr):
    # Clean up newlines
    expr = expr.replace("\n", "")

    # Convert -Kxxxxx to not Kxxxxx
    expr = re.sub(r'(?<!\w)-K(\d{5})', r'not K\1', expr)

    # Pad special chars with spaces
    expr = re.sub(r'([\(\)\+\-,])', r' \1 ', expr)
    expr = re.sub(r'\s+', ' ', expr).strip()

    # Replace KEGG-style operators with boolean logic
    expr = expr.replace('+', ' and ')
    expr = expr.replace(',', ' or ')
    expr = expr.replace('-', ' and not ')

    # Recursively fix adjacency "Kxxxxx Kyyyyy" => "Kxxxxx and Kyyyyy"
    while re.search(r'(K\d{5})\s+(K\d{5})', expr):
        expr = re.sub(r'(K\d{5})\s+(K\d{5})', r'\1 and \2', expr)

    expr = re.sub(r'(K\d{5})\s+\(', r'\1 and (', expr)
    expr = re.sub(r'\)\s+(K\d{5})', r') and \1', expr)
    expr = re.sub(r'\)\s+\(', r') and (', expr)

    return expr

def parse_and_filter(expr, valid_kos, mod_id):
    try:
        parsed = algebra.parse(expr)
        return recursively_filter(parsed, valid_kos)
    except Exception as e:
        print(f"❌ Failed to parse expression in {mod_id}: {expr}\nError: {e}")
        return None

def recursively_filter(node, valid_kos):
    if isinstance(node, boolean.Symbol):
        return node if node.obj in valid_kos else boolean.Symbol("KXXXXX")
    elif isinstance(node, boolean.NOT):
        return ~recursively_filter(node.args[0], valid_kos)
    elif isinstance(node, boolean.AND):
        return algebra.AND(*[recursively_filter(arg, valid_kos) for arg in node.args])
    elif isinstance(node, boolean.OR):
        return algebra.OR(*[recursively_filter(arg, valid_kos) for arg in node.args])
    return node

def stringify_filtered(expr):
    if expr == algebra.TRUE or expr == algebra.FALSE:
        return "KXXXXX"
    if isinstance(expr, boolean.Symbol):
        return expr.obj
    elif isinstance(expr, boolean.AND):
        parts = [stringify_filtered(arg) for arg in expr.args]
        return "(" + " and ".join(parts) + ")" if len(parts) > 1 else parts[0]
    elif isinstance(expr, boolean.OR):
        parts = [stringify_filtered(arg) for arg in expr.args]
        return "(" + " or ".join(parts) + ")" if len(parts) > 1 else parts[0]
    elif isinstance(expr, boolean.NOT):
        return f"(not {stringify_filtered(expr.args[0])})"
    return str(expr)

def extract_top_level_steps(def_string):
    """Extract top-level logical steps respecting parentheses, skipping '--'."""
    steps = []
    buffer = ''
    depth = 0
    i = 0

    while i < len(def_string):
        if def_string[i:i+2] == '--' and depth == 0:
            i += 2
            continue

        ch = def_string[i]
        if ch == '(':
            depth += 1
            buffer += ch
        elif ch == ')':
            depth -= 1
            buffer += ch
        elif ch == ' ' and depth == 0:
            if buffer.strip():
                steps.append(buffer.strip())
                buffer = ''
        else:
            buffer += ch
        i += 1

    if buffer.strip():
        steps.append(buffer.strip())

    return steps

def clean_and_validate_expression(expr_str):
    try:
        expr = algebra.parse(expr_str)
        simplified = expr.simplify()
        if simplified == algebra.FALSE:
            return None
        return stringify_filtered(simplified)
    except Exception:
        return None

def build_module_step_db():
    valid_kos = load_kos_from_hal(HAL_PATH)
    modules = extract_modules(KO_JSON)
    modules = sorted(modules, key=lambda x: int(x[0][1:]))

    with open(OUTPUT_MODULE_STEPS, "w") as out_f:
        out_f.write("name\tk_string\tmodule_step\n")
        for mod_id, name in tqdm(modules, desc="Modules", unit="module"):
            raw = get_kegg_module_definition(mod_id)
            if not raw:
                continue

            step_exprs = extract_top_level_steps(raw)
            # optional live postfix
            # tqdm.write(f"{mod_id} ({len(step_exprs)} steps)")  # less noisy than dynamic postfix

            for i, expr in enumerate(step_exprs, 1):
                normalized = normalize_expression(expr)
                parsed = parse_and_filter(normalized, valid_kos, mod_id)
                if parsed is None:
                    continue
                final_expr = stringify_filtered(parsed)
                cleaned_expr = clean_and_validate_expression(final_expr)
                if cleaned_expr is None:
                    continue
                out_f.write(f"{name}\t{cleaned_expr}\t{mod_id}+{i:02d}\n")

    print(f"✅ Wrote {OUTPUT_MODULE_STEPS}")

# -----------------------
# Main
# -----------------------
def main():
    build_hmm_table_template()
    collapse_hmm_table()
    build_module_step_db()

if __name__ == "__main__":
    main()
