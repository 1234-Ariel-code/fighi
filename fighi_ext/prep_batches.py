#!/usr/bin/env python3
import argparse, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top_cols", required=True, help="File: one SNP name per line (ranked)")
    ap.add_argument("--batch_size", type=int, default=1000)
    ap.add_argument("--halo", type=int, default=200, help="overlap between consecutive batches")
    ap.add_argument("--anchors", type=int, default=200, help="top-A anchors included in every batch")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cols = [ln.strip() for ln in open(args.top_cols) if ln.strip()]

    anchors = cols[:args.anchors]
    body = cols[args.anchors:]

    # step ensures: anchors + window + halo overlap across batches
    step = args.batch_size - args.halo - len(anchors)
    step = max(1, step)

    batches = []
    i = 0
    while i < len(body):
        window = body[i:i + (args.batch_size - len(anchors))]
        if not window:
            break
        batch = anchors + window
        batches.append(batch)
        i += step

    # write per-batch lists and master index
    lst = os.path.join(args.outdir, "batches.txt")
    with open(lst, "w") as f:
        for bi, batch in enumerate(batches, 1):
            path = os.path.join(args.outdir, f"batch_{bi:04d}.cols")
            with open(path, "w") as g:
                g.write("\n".join(batch) + "\n")
            f.write(path + "\n")
    print(f"[prep_batches] Wrote {len(batches)} batches in {args.outdir}")

if __name__ == "__main__":
    main()
