import re
from google.cloud.storage import Client, Bucket, Blob

# Regex for gs file locations: gs:://<bucket-name>/<blob-name>
GS_PATH_REGEX = re.compile(r'^gs://([a-z0-9._\-\s]+)((?:/[a-zA-Z0-9_.\-\s]+)*)$')


def get_bucket_and_blob_from_gs_path(storage_client: Client, gs_path: str) -> (Bucket, Blob):
    """
    Returns the bucket and the blob from a given google storage path.

    :param storage_client: the gcp storage client.
    :param gs_path: the path in the following format: (gs://<bucket-name>/<blob-name>).
    :return: a tuple containing both the bucket and blob: (bucket, blob).
    """
    match = GS_PATH_REGEX.match(gs_path)
    if not match:
        raise ValueError(f"'{gs_path}' is not a valid path!")

    groups = match.groups()
    bucket_name = groups[0]
    blob_name = groups[1][1:] if groups[1] else None  # the [1:] is to remove the leading slash from the blob name

    bucket: Bucket = storage_client.bucket(bucket_name=bucket_name)
    blob: Blob = bucket.get_blob(blob_name=blob_name)

    return bucket, blob


def get_file_name_from_blob_name(blob_name: str) -> str:
    """
    Given a blob name, it will treat it as a traditional file path and return the file name (the name after the last "/").

    :param blob_name: the full blob name, including (pseudo) subdirectories
    :return: the blobs file name
    """
    split = blob_name.split('/')
    return split[-1]
