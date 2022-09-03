import fnmatch
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

from google.cloud import storage

from gcp.base import GCPPath


@dataclass(frozen=True)
class GCPClient(object):
    client: storage.Client

    def file_exists(self, path: GCPPath) -> bool:
        return bool(self._get_blob(path).exists())

    def download_file(self, gcp_path: GCPPath, local_path: Path) -> None:
        logging.info(f"Starting download of '{gcp_path}' to '{local_path}'.")
        if not self.file_exists(gcp_path):
            raise FileNotFoundError(f"Cannot download file that doesn't exist: {gcp_path}")
        self._create_parent_dir_if_not_exists(local_path)
        self._get_blob(gcp_path).download_to_filename(str(local_path))
        if not local_path.exists():
            raise FileNotFoundError(f"Download of '{gcp_path}' to '{local_path}' has failed.")
        logging.info(f"Finished download of '{gcp_path}' to '{local_path}'.")

    def upload_file(self, local_path: Path, gcp_path: GCPPath) -> None:
        logging.info(f"Starting upload of '{local_path}' to '{gcp_path}'.")
        if not local_path.exists():
            raise FileNotFoundError(f"Cannot upload file that doesn't exist: '{local_path}'")
        self._get_blob(gcp_path).upload_from_filename(str(local_path))
        if not self.file_exists(gcp_path):
            raise FileNotFoundError(f"Upload of '{local_path}' to '{gcp_path}' has failed.")
        logging.info(f"Finished upload of '{local_path}' to '{gcp_path}'.")

    def get_files_in_directory(self, path: GCPPath) -> List[GCPPath]:
        if path.relative_path[-1] != "/":
            prefix = path.relative_path + "/"
        else:
            prefix = path.relative_path
        blobs = self.client.list_blobs(path.bucket_name, prefix=prefix, delimiter="/")
        return [GCPPath(path.bucket_name, blob.name) for blob in blobs]

    def get_matching_file_paths(self, path: GCPPath) -> List[GCPPath]:
        matching_paths: List[GCPPath] = []
        prefix_to_match = path.relative_path.split("*")[0]
        for blob in self.client.list_blobs(path.bucket_name, prefix=prefix_to_match):
            if fnmatch.fnmatch(blob.name, path.relative_path):
                matching_paths.append(GCPPath(path.bucket_name, blob.name))
        return matching_paths

    def get_text(self, path: GCPPath) -> str:
        return self._get_blob(path).download_as_text()

    def _get_blob(self, path: GCPPath) -> storage.Blob:
        return self.client.bucket(path.bucket_name).blob(path.relative_path)

    def _create_parent_dir_if_not_exists(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)