from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class GCPPath(object):
    bucket_name: str
    relative_path: str

    @classmethod
    def from_string(cls, path: str) -> "GCPPath":
        if not path.startswith("gs://"):
            raise ValueError(f"Path is not a GCP bucket path: '{path}'")
        bucket_name = path.split("/")[2]
        relative_path = "/".join(path.split("/")[3:])
        return GCPPath(bucket_name, relative_path)

    def __str__(self) -> str:
        return f"gs://{self.bucket_name}/{self.relative_path}"

    def get_parent_directory(self) -> "GCPPath":
        if self.relative_path[-1] == "/":
            relative_path = self.relative_path[:-1]
        else:
            relative_path = self.relative_path
        parent_directory_relative_path = "/".join(relative_path.split("/")[:-1])
        return GCPPath(self.bucket_name, parent_directory_relative_path)

    def append_suffix(self, suffix: str) -> "GCPPath":
        return GCPPath(self.bucket_name, self.relative_path + suffix)
