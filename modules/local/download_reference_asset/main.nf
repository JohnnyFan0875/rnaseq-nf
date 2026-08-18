process download_reference_asset {

    tag "${asset_name}"

    container 'python:3.11-slim'

    input:
    val asset_name
    val asset_url
    val auto_download

    output:
    path "${asset_name}", emit: asset
    path "${asset_name}", emit: gtf, optional: true
    path "${asset_name}", emit: fasta, optional: true

    script:
    """
    if [ "${auto_download}" != "true" ]; then
        echo "ERROR: Missing reference asset ${asset_name} and reference_auto_download=false." >&2
        exit 1
    fi

    python - <<'PY'
import gzip
import shutil
import urllib.request

url = "${asset_url}"
out_path = "${asset_name}"
tmp_path = out_path + ".gz"

with urllib.request.urlopen(url) as response, open(tmp_path, "wb") as out_handle:
    shutil.copyfileobj(response, out_handle)

with gzip.open(tmp_path, "rb") as src, open(out_path, "wb") as dst:
    shutil.copyfileobj(src, dst)
PY
    rm -f "${asset_name}.gz"
    """
}
