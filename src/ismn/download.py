import requests
import zipfile
import zlib
from pathlib import Path
from ismn.misc import collect_stm_cov, write_overview
from ismn.const import nrt_networks

class ISMNDownloader:
    BASE_URL = "https://ismn.earth"
    LOGIN_URL = f"{BASE_URL}/en/accounts/login/"
    DOWNLOAD_URL = f"{BASE_URL}/en/dataviewer/api/download_archive"

    def __init__(self, username: str, password: str, output_path: str = "ismn_archive.zip"):
        self.username = username
        self.password = password
        self.output_path = Path(output_path).expanduser().resolve()
        self.session = requests.Session()

    def _get_csrf_token(self) -> str:
        print("Step 1: Fetching login page and CSRF token...")
        self.session.get(self.LOGIN_URL)
        csrf_token = self.session.cookies["csrftoken"]
        print(f"  ✓ CSRF token obtained: {csrf_token[:10]}...")
        return csrf_token

    def _login(self, csrf_token: str) -> None:
        print("\nStep 2: Logging in...")
        response = self.session.post(
            self.LOGIN_URL,
            data={
                "csrfmiddlewaretoken": csrf_token,
                "login": self.username,
                "password": self.password,
            },
            headers={"Referer": self.LOGIN_URL},
        )
        if response.ok:
            print(f"  ✓ Login successful (status {response.status_code})")
        else:
            raise RuntimeError(f"Login failed (status {response.status_code})")

    def _download(self) -> None:
        print("\nStep 3: Starting archive download...")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        response = self.session.get(self.DOWNLOAD_URL, stream=True)

        if not response.ok:
            raise RuntimeError(f"Download request failed (status {response.status_code})")

        total_size = int(response.headers.get("content-length", 0))
        if total_size:
            print(f"  ✓ File size: {total_size / (1024**3):.2f} GB")
        else:
            print("  ! File size unknown (no Content-Length header)")

        chunk_size = 1024 * 1024  # 1 MB
        downloaded = 0

        with open(self.output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size:
                        percent = downloaded / total_size * 100
                        downloaded_gb = downloaded / (1024**3)
                        total_gb = total_size / (1024**3)
                        print(f"  Downloading... {downloaded_gb:.2f} / {total_gb:.2f} GB ({percent:.1f}%)", end="\r")
                    else:
                        print(f"  Downloaded: {downloaded / (1024**3):.2f} GB", end="\r")

        print(f"\n  ✓ Download complete: {self.output_path} ({downloaded / (1024**3):.2f} GB)")

    def run(self) -> None:
        csrf_token = self._get_csrf_token()
        self._login(csrf_token)
        self._download()



class ISMNExtractor:
    def __init__(self, archive_path: str, nrt_networks: list[str], output_dir: str = "ISMN"):
        self.archive_path = Path(archive_path)
        self.nrt_networks = nrt_networks
        self.output_dir = Path(output_dir)

    @staticmethod
    def _file_crc32(path: Path, chunk_size: int = 65536) -> int:
        crc = 0
        with open(path, "rb") as f:
            while chunk := f.read(chunk_size):
                crc = zlib.crc32(chunk, crc)
        return crc

    def _needs_extraction(self, zinfo: zipfile.ZipInfo, target: Path) -> tuple[bool, str]:
        """Return (needs_extraction, reason)."""
        if not target.exists():
            return True, "new"
        if target.stat().st_size != zinfo.file_size:
            return True, "size differs"
        if self._file_crc32(target) != zinfo.CRC:
            return True, "content differs"
        return False, "unchanged"

    def run(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting {len(self.nrt_networks)} NRT networks from {self.archive_path}...")
        print(f"Output directory: {self.output_dir.resolve()}\n")

        totals = {"new": 0, "updated": 0, "skipped": 0, "missing_networks": 0}

        with zipfile.ZipFile(self.archive_path, "r") as zf:
            all_entries = zf.namelist()

            for network in self.nrt_networks:
                matching = [e for e in all_entries if e.startswith(f"{network}/")]

                if not matching:
                    print(f"  ! Network not found in archive: {network}")
                    totals["missing_networks"] += 1
                    continue

                print(f"\n  Processing {network} ({len(matching)} files)...")
                net_new = net_updated = net_skipped = 0

                for entry in matching:
                    # Directory entries in zips end with "/"; skip them
                    if entry.endswith("/"):
                        continue

                    zinfo = zf.getinfo(entry)
                    target = self.output_dir / entry
                    needs, reason = self._needs_extraction(zinfo, target)

                    if needs:
                        zf.extract(entry, self.output_dir)
                        if reason == "new":
                            net_new += 1
                            print(f"    + new      {entry}")
                        else:
                            net_updated += 1
                            print(f"    ~ updated  {entry}  ({reason})")
                    else:
                        net_skipped += 1
                        print(f"    = skipped  {entry}")

                totals["new"] += net_new
                totals["updated"] += net_updated
                totals["skipped"] += net_skipped
                print(f"  ✓ {network}: {net_new} new, {net_updated} updated, {net_skipped} unchanged")

        print(
            f"\n✓ All done. "
            f"{totals['new']} new, {totals['updated']} updated, {totals['skipped']} unchanged."
        )
        if totals["missing_networks"]:
            print(f"  ! {totals['missing_networks']} requested network(s) not found in archive.")


        print("Collecting coverage information from .stm files...")
        period_to = collect_stm_cov(str(self.output_dir), n_proc=4) 
        print("Writing overview file...")  
        write_overview(self.output_dir, period_to=period_to)


if __name__ == "__main__":
    ISMNExtractor(archive_path="ismn_archive.zip", nrt_networks=nrt_networks, output_dir="ISMN").run()