#!/usr/bin/env python3
"""Build the containers of a local module with Wave (Seqera Containers).

Sends the module's environment.yml to the public Wave service twice (Docker
image and Singularity image), waits for both builds and prints the container
directive. With --write, the directive replaces the one in the module main.nf.

The images are frozen in community.wave.seqera.io: public, permanent and
identified by a hash of the environment, so an unchanged environment gives the
same images. No account or token is needed.

Usage: dev/nfdev container modules/local/<tool>/<subtool> [--write]
"""

import base64
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

WAVE = "https://wave.seqera.io"


def request(url, payload=None):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)


def submit(environment, image_format):
    payload = {
        "packages": {
            "type": "CONDA",
            "environment": base64.b64encode(environment.encode()).decode(),
            "channels": ["conda-forge", "bioconda"],
        },
        "freeze": True,
        "format": image_format,
        "containerPlatform": "linux/amd64",
    }
    return request(f"{WAVE}/v1alpha2/container", payload)


def wait(build_id):
    while True:
        status = request(f"{WAVE}/v1alpha1/builds/{build_id}/status")
        if status.get("status") == "COMPLETED":
            return status.get("succeeded", False)
        time.sleep(20)


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    module = Path(sys.argv[1])
    write = "--write" in sys.argv[2:]
    environment = (module / "environment.yml").read_text()

    images = {}
    for image_format in ("docker", "sif"):
        build = submit(environment, image_format)
        images[image_format] = build["containerImage"]
        if build.get("cached") and build.get("succeeded") is False:
            sys.exit(f"{image_format}: an earlier build of this environment failed (cached): {build.get('buildId')}")
        print(f"{image_format}: {build['containerImage']} (build {build.get('buildId')}, cached: {build.get('cached')})", flush=True)
        if build.get("buildId") and not build.get("cached"):
            print(f"{image_format}: building, see {WAVE}/view/builds/{build['buildId']}", flush=True)
            if not wait(build["buildId"]):
                sys.exit(f"{image_format}: build failed, see {WAVE}/view/builds/{build['buildId']}")

    directive = (
        "    container \"${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?\n"
        f"        '{images['sif']}' :\n"
        f"        '{images['docker']}' }}\"\n"
    )
    print("\n" + directive)
    if write:
        main_nf = module / "main.nf"
        text = main_nf.read_text()
        pattern = re.compile(r"    container (?:'[^']*'|\"\$\{[^}]*\}\")\n", re.S)
        if not pattern.search(text):
            sys.exit(f"no container directive found in {main_nf}")
        main_nf.write_text(pattern.sub(lambda _: directive, text, count=1))
        print(f"updated {main_nf}")


if __name__ == "__main__":
    main()
