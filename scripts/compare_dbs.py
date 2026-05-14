#!/usr/bin/env python3
"""Database utilities: compare IDs across two databases, or reindex a database to match a reference."""

import argparse
import os
import sys


def load_index(index_path):
    """Load index file into a dict: id -> (start, length)."""
    lookup = {}
    with open(index_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            entry_id, start, length = parts[0], int(parts[1]), int(parts[2])
            lookup[entry_id] = (start, length)
    return lookup


def read_record(db_fh, start, length):
    """Read a single record from the data file."""
    db_fh.seek(start)
    return db_fh.read(length)


def compare_ids(args):
    """For each ID in db1, check if db2 has the same data for that ID."""
    print(f"Loading {args.db2} index ...")
    db2_index = load_index(args.db2 + ".index")
    print(f"  {len(db2_index)} entries")

    print(f"Streaming {args.db1} index and comparing ...")
    total = 0
    shared = 0
    mismatches = 0
    example_shown = False

    with open(args.db1, "rb") as db1_fh, open(args.db2, "rb") as db2_fh:
        with open(args.db1 + ".index") as idx:
            for line in idx:
                parts = line.strip().split("\t")
                entry_id, start1, length1 = parts[0], int(parts[1]), int(parts[2])
                total += 1

                if entry_id not in db2_index:
                    continue

                shared += 1
                start2, length2 = db2_index[entry_id]
                rec1 = read_record(db1_fh, start1, length1)
                rec2 = read_record(db2_fh, start2, length2)

                if rec1 != rec2:
                    mismatches += 1
                    if not example_shown:
                        print(f"\nExample mismatch for ID '{entry_id}':")
                        print(f"  db1: {rec1!r}")
                        print(f"  db2: {rec2!r}")
                        example_shown = True

    print(f"\nTotal IDs in db1: {total}")
    print(f"Shared IDs: {shared}")
    print(f"Mismatches: {mismatches}")
    if not mismatches:
        print("All shared IDs have identical records.")


def reindex(args):
    """Reindex db to use the same numeric IDs as ref_db, matched by header content."""
    ref_db = args.ref_db
    db = args.db
    db_out = args.db_out

    # Step 1: Build reverse index for db: header_content -> db_numeric_id
    print(f"Building reverse index from {db}_h ...")
    db_h_index = load_index(db + "_h.index")
    reverse = {}  # header bytes -> db numeric id
    with open(db + "_h", "rb") as db_h_fh:
        for entry_id, (start, length) in db_h_index.items():
            header = read_record(db_h_fh, start, length)
            if header in reverse:
                print(f"Error: duplicate header in {db}_h for IDs {reverse[header]} and {entry_id}", file=sys.stderr)
                sys.exit(1)
            reverse[header] = entry_id
    print(f"  {len(reverse)} entries")

    # Step 2: Load db indices
    db_index = load_index(db + ".index")
    db_h_index_reload = load_index(db + "_h.index")
    has_ss = os.path.exists(db + "_ss") and os.path.exists(db + "_ss.index")
    if has_ss:
        db_ss_index = load_index(db + "_ss.index")
        print(f"  Found {db}_ss ({len(db_ss_index)} entries)")
    has_ca = os.path.exists(db + "_ca") and os.path.exists(db + "_ca.index")
    if has_ca:
        db_ca_index = load_index(db + "_ca.index")
        print(f"  Found {db}_ca ({len(db_ca_index)} entries)")

    # Step 3: Stream ref_db_h index, look up each header in reverse index, write new indices
    print(f"Loading {ref_db} index ...")
    ref_db_index = load_index(ref_db + ".index")
    print(f"  {len(ref_db_index)} entries")

    print(f"Mapping {ref_db} IDs to {db} ...")
    out_index_lines = []
    out_h_index_lines = []
    out_ss_index_lines = []
    out_ca_index_lines = []
    skipped = 0

    with open(ref_db + "_h", "rb") as ref_h_fh:
        with open(ref_db + "_h.index") as ref_h_idx:
            for line in ref_h_idx:
                parts = line.strip().split("\t")
                ref_id, start, length = parts[0], int(parts[1]), int(parts[2])

                if ref_id not in ref_db_index:
                    skipped += 1
                    continue

                header = read_record(ref_h_fh, start, length)

                if header not in reverse:
                    header_text = header.rstrip(b"\n\0").decode(errors="replace")
                    print(f"Error: ref_db entry {ref_id} (header '{header_text}') not found in {db}", file=sys.stderr)
                    sys.exit(1)

                db_id = reverse[header]

                # New index entry: ref_db's numeric ID -> db's data offset
                db_start, db_length = db_index[db_id]
                out_index_lines.append(f"{ref_id}\t{db_start}\t{db_length}\n")

                # Same for _h
                db_h_start, db_h_length = db_h_index_reload[db_id]
                out_h_index_lines.append(f"{ref_id}\t{db_h_start}\t{db_h_length}\n")

                # Same for _ss if present
                if has_ss:
                    db_ss_start, db_ss_length = db_ss_index[db_id]
                    out_ss_index_lines.append(f"{ref_id}\t{db_ss_start}\t{db_ss_length}\n")

                # Same for _ca if present
                if has_ca:
                    db_ca_start, db_ca_length = db_ca_index[db_id]
                    out_ca_index_lines.append(f"{ref_id}\t{db_ca_start}\t{db_ca_length}\n")

    # Step 4: Create symlinks and write new index files
    db = os.path.abspath(db)
    db_out = os.path.abspath(db_out)

    def force_symlink(target, link_name):
        if os.path.islink(link_name) or os.path.exists(link_name):
            os.remove(link_name)
        os.symlink(target, link_name)

    force_symlink(db, db_out)
    force_symlink(db + ".dbtype", db_out + ".dbtype")
    force_symlink(db + "_h", db_out + "_h")
    force_symlink(db + "_h.dbtype", db_out + "_h.dbtype")

    with open(db_out + ".index", "w") as f:
        f.writelines(out_index_lines)
    with open(db_out + "_h.index", "w") as f:
        f.writelines(out_h_index_lines)

    if has_ss:
        force_symlink(db + "_ss", db_out + "_ss")
        force_symlink(db + "_ss.dbtype", db_out + "_ss.dbtype")
        with open(db_out + "_ss.index", "w") as f:
            f.writelines(out_ss_index_lines)

    if has_ca:
        force_symlink(db + "_ca", db_out + "_ca")
        force_symlink(db + "_ca.dbtype", db_out + "_ca.dbtype")
        with open(db_out + "_ca.index", "w") as f:
            f.writelines(out_ca_index_lines)

    print(f"  Wrote {len(out_index_lines)} entries (skipped {skipped} ref_db_h-only entries)")
    print(f"  {db_out} -> {db}")
    print(f"  {db_out}_h -> {db}_h")
    if has_ss:
        print(f"  {db_out}_ss -> {db}_ss")
    if has_ca:
        print(f"  {db_out}_ca -> {db}_ca")
    print("Done.")


def load_lookup(path):
    """Load a .lookup file (numeric_id<TAB>name<TAB>file_idx) into name -> numeric_id."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out[parts[1]] = int(parts[0])
    return out


def build_prefdb(args):
    """Build a single-shard prefDB from a 2-col tsv (qname<TAB>tname).
    Only the target-id column matters to structurealign; score/diag are placeholders.
    """
    qmap = load_lookup(args.qdb + ".lookup")

    # Collect needed target names from tsv, plus first-hit per query.
    first_hit = {}  # qname -> tname (first occurrence only)
    need_t = set()
    with open(args.tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            q, t = parts[0], parts[1]
            if q in first_hit:
                continue
            first_hit[q] = t
            need_t.add(t)

    # Stream tdb.lookup for only the names we need.
    tmap = {}
    with open(args.tdb + ".lookup") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if parts[1] in need_t:
                tmap[parts[1]] = int(parts[0])
                if len(tmap) == len(need_t):
                    break
    missing_t = need_t - set(tmap)
    if missing_t:
        ex = next(iter(missing_t))
        print(f"Error: {len(missing_t)} target names not in {args.tdb}.lookup (e.g. {ex!r})", file=sys.stderr)
        sys.exit(1)

    # Build (qid, record_bytes) sorted by qid, write data + index.
    # Query names not in qdb are silently skipped: e.g. when the search ran
    # over the full pred db but qdb is a subset for quicker testing.
    records = []
    skipped_q = 0
    for qname, tname in first_hit.items():
        if qname not in qmap:
            skipped_q += 1
            continue
        rec = f"{tmap[tname]}\t255\t0\n\0".encode()
        records.append((qmap[qname], rec))
    records.sort(key=lambda x: x[0])

    with open(args.out, "wb") as f_data, open(args.out + ".index", "w") as f_idx:
        offset = 0
        for qid, rec in records:
            f_data.write(rec)
            f_idx.write(f"{qid}\t{offset}\t{len(rec)}\n")
            offset += len(rec)
    with open(args.out + ".dbtype", "wb") as f:
        f.write(b"\x07\x00\x00\x00")  # PREFILTER_RES
    print(f"Wrote prefDB with {len(records)} queries to {args.out} (skipped {skipped_q} queries not in qdb)")


def main():
    parser = argparse.ArgumentParser(description="Database utilities")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_compare = subparsers.add_parser("compare_ids", help="Compare records for shared IDs across two databases")
    p_compare.add_argument("db1", help="Path to first database file (without .index)")
    p_compare.add_argument("db2", help="Path to second database file (without .index)")
    p_compare.set_defaults(func=compare_ids)

    p_reindex = subparsers.add_parser("reindex", help="Reindex db to match ref_db numbering via header matching")
    p_reindex.add_argument("ref_db", help="Reference database (provides the target numeric IDs)")
    p_reindex.add_argument("db", help="Source database (data to reindex)")
    p_reindex.add_argument("db_out", help="Output database path (symlinks to db, new index)")
    p_reindex.set_defaults(func=reindex)

    p_build = subparsers.add_parser("build_prefdb", help="Build a prefDB from a 2-col tsv (qname<TAB>tname)")
    p_build.add_argument("tsv", help="Input tsv with qname<TAB>tname; first row per qname is used")
    p_build.add_argument("qdb", help="Query db (its .lookup maps qname -> numeric id)")
    p_build.add_argument("tdb", help="Target db (its .lookup maps tname -> numeric id)")
    p_build.add_argument("out", help="Output prefDB path (writes <out>, <out>.index, <out>.dbtype)")
    p_build.set_defaults(func=build_prefdb)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
