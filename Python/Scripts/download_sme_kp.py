#!/usr/bin/env python3
"""Download the Kp index (GFZ Potsdam) and the SME/SML/SMU indices (SuperMAG)
for a given UTC time period and save each to its own CSV file.

Only the Python standard library is used, so no extra packages are needed.

Kp is fetched from GFZ Potsdam's public web service. No registration is
required.

SME/SML/SMU are SuperMAG-exclusive derived indices, so they can only be
obtained from SuperMAG's web service. That service requires a free "logon"
username: register once (just an email address, no API key) at
https://supermag.jhuapl.edu, then pass the username with --userid or set it
in the SUPERMAG_USERID environment variable.

Examples:
    python download_sme_kp.py --start 2024-05-10T00:00:00 \\
        --end 2024-05-13T00:00:00 --userid myname --out period1

    SUPERMAG_USERID=myname python download_sme_kp.py \\
        --start 2024-05-10T00:00:00 --end 2024-05-13T00:00:00 --out period1

Run once per time period; it writes:
    <out>_kp.csv   3-hour Kp index
    <out>_sme.csv  1-minute SME, SML, SMU index
"""

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

GFZ_KP_URL = "https://kp.gfz.de/app/json/"
SUPERMAG_INDICES_URL = "https://supermag.jhuapl.edu/services/indices.php"
TIMEOUT_SEC = 60


def parse_utc(value: str) -> datetime:
    """Parse an ISO-8601 timestamp (optionally with a trailing Z) as UTC."""
    return datetime.fromisoformat(value.replace("Z", "+00:00")).replace(
        tzinfo=timezone.utc
    )


def fetch_json(url: str) -> dict | list:
    req = Request(url, headers={"User-Agent": "download_sme_kp.py"})
    try:
        with urlopen(req, timeout=TIMEOUT_SEC) as resp:
            raw_text = resp.read().decode("utf-8", errors="replace")
    except (HTTPError, URLError) as exc:
        raise RuntimeError(f"Could not reach {url}: {exc}") from exc

    if raw_text.startswith("ERROR"):
        raise RuntimeError(f"Server returned an error for {url}:\n{raw_text}")
    return json.loads(raw_text)


def download_kp(start: datetime, end: datetime, out_csv: str) -> None:
    """Fetch the 3-hour Kp index from GFZ Potsdam (no registration needed)."""
    params = urlencode({
        "start": start.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "end": end.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "index": "Kp",
        "status": "def",
    })
    data = fetch_json(f"{GFZ_KP_URL}?{params}")

    times = data["datetime"]
    kp_values = data["Kp"]
    statuses = data.get("status", [None] * len(times))

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["datetime_utc", "Kp", "status"])
        for dt_str, kp_val, status in zip(times, kp_values, statuses):
            writer.writerow([dt_str, kp_val, status])

    print(f"  Saved {len(times)} Kp records to {out_csv}")


def download_sme(userid: str, start: datetime, end: datetime, out_csv: str) -> None:
    """Fetch 1-minute SME/SML/SMU from SuperMAG (requires a free logon)."""
    extent_sec = int((end - start).total_seconds())
    if extent_sec <= 0:
        raise ValueError("--end must be after --start")

    params = urlencode({
        "python": "",
        "nohead": "",
        "start": start.strftime("%Y-%m-%dT%H:%M"),
        "logon": userid,
        "extent": "%12.12d" % extent_sec,
        "indices": "sme,sml,smu",
    })
    # urlencode would otherwise turn the bare "python"/"nohead" flags into
    # "python=" which the SuperMAG service does not accept.
    url = f"{SUPERMAG_INDICES_URL}?{params}".replace("python=&", "python&").replace(
        "nohead=&", "nohead&"
    )
    records = fetch_json(url)

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["datetime_utc", "SME", "SML", "SMU"])
        for row in records:
            dt_str = datetime.fromtimestamp(
                row["tval"], tz=timezone.utc
            ).strftime("%Y-%m-%dT%H:%M:%SZ")
            writer.writerow([dt_str, row.get("SME"), row.get("SML"), row.get("SMU")])

    print(f"  Saved {len(records)} SME/SML/SMU records to {out_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--start", required=True, help="start time, UTC ISO format, e.g. 2024-05-10T00:00:00"
    )
    parser.add_argument(
        "--end", required=True, help="end time, UTC ISO format, e.g. 2024-05-13T00:00:00"
    )
    parser.add_argument(
        "--out", required=True, help="output file prefix, e.g. 'period1'"
    )
    parser.add_argument(
        "--userid",
        default=os.environ.get("SUPERMAG_USERID"),
        help="SuperMAG logon username (or set SUPERMAG_USERID env var); "
        "register for free at https://supermag.jhuapl.edu. "
        "If omitted, SME/SML/SMU are skipped.",
    )
    args = parser.parse_args()

    start = parse_utc(args.start)
    end = parse_utc(args.end)

    print(f"Downloading Kp for {args.start} to {args.end}...")
    download_kp(start, end, f"{args.out}_kp.csv")

    if args.userid:
        print(f"Downloading SME/SML/SMU for {args.start} to {args.end}...")
        download_sme(args.userid, start, end, f"{args.out}_sme.csv")
    else:
        print(
            "No --userid / SUPERMAG_USERID given, skipping SME/SML/SMU download.",
            file=sys.stderr,
        )

    print("Done!")


if __name__ == "__main__":
    main()
