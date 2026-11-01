import modal

app = modal.App("her2-scfv-folding")

image = modal.Image.debian_slim(python_version="3.11").pip_install("boltz")
cache = modal.Volume.from_name("boltz-cache", create_if_missing=True)


@app.function(image=image, gpu="L40S", timeout=3600, volumes={"/root/.boltz": cache})
def predict(yaml_text: str, name: str):
    import pathlib
    import subprocess

    work = pathlib.Path("/tmp/work")
    work.mkdir(exist_ok=True)
    (work / f"{name}.yaml").write_text(yaml_text)
    result = subprocess.run(
        ["boltz", "predict", str(work / f"{name}.yaml"),
         "--out_dir", str(work / "out"), "--use_msa_server",
         "--write_full_pae"],
        capture_output=True, text=True,
    )
    # Always show Boltz's own messages so that errors are visible
    print("=== BOLTZ STDOUT (end) ===")
    print(result.stdout[-4000:])
    print("=== BOLTZ STDERR (end) ===")
    print(result.stderr[-8000:])
    if result.returncode != 0:
        raise RuntimeError(f"boltz failed with exit code {result.returncode}")
    cache.commit()
    files = {}
    for p in (work / "out").rglob("*"):
        if p.suffix in (".cif", ".json", ".npz"):
            files[str(p.relative_to(work / "out"))] = p.read_bytes()
    return files


@app.local_entrypoint()
def main(yaml_path: str):
    import pathlib

    src = pathlib.Path(yaml_path)
    results = predict.remote(src.read_text(), src.stem)
    out = pathlib.Path("results") / src.stem
    for rel, data in results.items():
        target = out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        print("saved", target)
