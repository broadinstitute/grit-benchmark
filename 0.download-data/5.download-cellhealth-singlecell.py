"""
Downloading single cell Cell Painting profiles from the Cell Health experiment

Note, together the files are 130 GB, and downloading will take a long time.

See https://github.com/broadinstitute/cell-health/blob/master/0.download-data/download.py
"""

import pathlib
import requests


def download_sqllite_file(filename, url):
    with requests.get(url, stream=True) as sql_request:
        sql_request.raise_for_status()
        with open(filename, "wb") as sql_fh:
            for chunk in sql_request.iter_content(chunk_size=819200000):
                if chunk:
                    assert isinstance(chunk, object)
                    sql_fh.write(chunk)


file_info = {
    "SQ00014610": "18028784",
    "SQ00014611": "18508583",
    "SQ00014612": "18505937",
    "SQ00014613": "18506036",
    "SQ00014614": "18031619",
    "SQ00014615": "18506108",
    "SQ00014616": "18506912",
    "SQ00014617": "18508316",
    "SQ00014618": "18508421",
}

download_dir = pathlib.Path("data/cell_health")
download_dir.mkdir(exist_ok=True, parents=True)

for plate in file_info:
    figshare_id = file_info[plate]
    filename = pathlib.Path(download_dir, f"{plate}.sqlite")
    if filename.exists():
        continue
    print(f"Now downloading... {filename}")
    url = f"https://nih.figshare.com/ndownloader/files/{figshare_id}"
    download_sqllite_file(filename, url)
    print("Done...\n\n")
