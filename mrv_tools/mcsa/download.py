import requests
import math
import multiprocessing as mp
from bs4 import BeautifulSoup
import os
from mrv_tools.tools import message, warning_m


def request_entry_data(verbose: int, njobs: int):
    url = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json'
    message(f"Requesting base page entries from: {url}", verbose=verbose)
    req = requests.get(url).json()
    num_entries = req["count"]
    num_entries_per_page = len(req["results"])
    num_pages = math.ceil(num_entries/num_entries_per_page)
    page_urls = [f"{url}&page={i+2}" for i in range(num_pages-1)]

    message(f"Total entries: {num_entries}, entries per page: {num_entries_per_page}, pages: {num_pages}",
            verbose=verbose)

    with mp.Pool(processes=njobs) as pool:
        page_requests = pool.map(requests.get, page_urls)

    results = req["results"]
    for page_req in page_requests:
        results.extend(page_req.json()["results"])
    assert(num_entries == len(results))
    return results


def _prepare_entry_arguments(_entries, root_dir, verbose: int = 0):
    urls = []
    for e in _entries:
        for m in e["reaction"]["mechanisms"]:
            if m["is_detailed"] is False:
                continue
            for s in m["steps"]:
                if s["is_product"] is True:
                    continue
                mcsa_full_id = f"{e['mcsa_id']}-{m['mechanism_id']}-{s['step_id']}"
                xml_url = s['marvin_xml']
                if len(xml_url) == 0:
                    warning_m(f"No URL found for \"{mcsa_full_id}\"", verbose=verbose)
                    continue
                urls.append(
                    (mcsa_full_id, xml_url, root_dir)
                )
    return urls


def _download_mrv_file(step_id, url, root_dir):
    if not url.startswith("http://"):
        url = "http://" + url

    implicit_path = os.path.join(root_dir, f"{step_id}.mrv")

    # Download xml data
    content = requests.get(url).content

    # make the files prettier
    prettier = BeautifulSoup(content, "xml").prettify()
    with open(implicit_path, "w") as f:
        f.write(prettier)

    return True


def request_marvin_files(root_dir: str, _entries, verbose: int, njobs: int):
    message(f"Requesting marvin files for {len(_entries)} M-CSA entries.", verbose=verbose)

    root_dir = os.path.join(root_dir, "implicit")
    if not os.path.exists(root_dir):
        os.makedirs(root_dir)

    mrv_args = _prepare_entry_arguments(_entries, root_dir, verbose=verbose)

    message(f"Found {len(mrv_args)} marvin file URLs", verbose=verbose)
    message("Downloading marvin files.", verbose=verbose)
    with mp.Pool(processes=njobs) as pool:
        pool.starmap(_download_mrv_file, mrv_args)


if __name__ == "__main__":
    entries = request_entry_data(True, 10)
    message(f"Found {len(entries)}")
    request_marvin_files("test", entries, True, 4)
