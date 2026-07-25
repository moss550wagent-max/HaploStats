#!/usr/bin/env python3
"""
HaploStats API Test Script
Sends the Phase 4 ambiguous mock patient to the running FastAPI server.
"""

import json
import sys
import time

try:
    import requests
except ImportError:
    print("❌ requests not installed. Run: pip install requests")
    sys.exit(1)

BASE_URL = "http://localhost:8000"

# Phase 4 mock patient: typed only at A, B, DRB1
MOCK_PATIENT = {
    "hla_a":      ["A*02:01:01:01", "A*01:01:01:01"],
    "hla_c":      None,
    "hla_b":      ["B*44:02:01:01", "B*08:01:01:01"],
    "hla_drb345": None,
    "hla_drb1":   ["DRB1*04:01:01:01SG", "DRB1*03:01:01:01SG"],
    "hla_dqa1":   None,
    "hla_dqb1":   None,
    "hla_dpa1":   None,
    "hla_dpb1":   None,
}


def test_health():
    print("\n🏥 GET /health")
    resp = requests.get(f"{BASE_URL}/health", timeout=10)
    print(f"   Status: {resp.status_code}")
    data = resp.json()
    print(f"   Response: {json.dumps(data, indent=2)}")
    assert resp.status_code == 200
    assert data["status"] == "ok"
    return True


def test_impute():
    print("\n📤 POST /impute  (Phase 4 ambiguous mock patient)")
    print(f"   Payload: {json.dumps(MOCK_PATIENT, indent=2)}")

    resp = requests.post(
        f"{BASE_URL}/impute?population=Global",
        json=MOCK_PATIENT,
        timeout=60,
        headers={"Content-Type": "application/json"},
    )

    print(f"\n   Status: {resp.status_code}")

    if resp.status_code != 200:
        print(f"   ❌ Error: {resp.text}")
        return False

    data = resp.json()

    print(f"\n{'=' * 72}")
    print(f"📦 RAW API RESPONSE")
    print(f"{'=' * 72}")
    print(json.dumps(data, indent=2))

    # Basic assertions
    assert data["status"] == "success"
    assert data["total_possible_pairs"] > 0
    assert len(data["top_pairs"]) == 3
    assert len(data["blocks"]) == 3

    print(f"\n{'=' * 72}")
    print(f"✅ All assertions passed!")
    print(f"{'=' * 72}")

    return True


def main():
    print("=" * 72)
    print("  HaploStats — API Test Suite")
    print("=" * 72)

    # Wait for server to be ready
    for attempt in range(10):
        try:
            r = requests.get(f"{BASE_URL}/health", timeout=3)
            if r.status_code == 200:
                print("  ✅ Server is up.\n")
                break
        except requests.exceptions.ConnectionError:
            if attempt == 0:
                print("  ⏳ Waiting for server to start...")
            time.sleep(1)
    else:
        print("  ❌ Could not connect to server. Is uvicorn running?")
        print(f"     Start it: uvicorn scripts.api:app --host 0.0.0.0 --port 8000")
        sys.exit(1)

    success = True
    success &= test_health()
    success &= test_impute()

    if success:
        print("\n🎉 All tests passed!")
    else:
        print("\n❌ Some tests failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
