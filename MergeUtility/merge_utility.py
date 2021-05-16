import os
import shutil
import argparse
from pathlib import Path
from subprocess import Popen
from collections import defaultdict


class MergeLanes:
    def __init__(self, in_dir: Path, out_dir: Path):
        self.in_dir = in_dir
        self.out_dir = out_dir

    # Group files together based on sample-name
    def group_fastqs(self, in_dir: Path, out_dir: Path) -> dict:
        """Function groups fastq files
        together to return a file group
        Args:
            in_dir - Input folder with bclconvert
                        output
            out_dir - Output folder[/FASTQ]
        """
        ext = "001.fastq.gz"  # default ext for fastq files
        files = sorted(os.listdir(in_dir))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        file_groups: dict = defaultdict(lambda: [])
        for file in files:
            if file.endswith(ext):
                sam_name = "_".join(file.split("_")[:-4])
                sam_num = file.split("_")[-4]
                sam_stn = file.split("_")[-2]
                key = f"{sam_name}_{sam_num}_{sam_stn}_{ext}"
                file_groups[key].append(f"{in_dir}/{file}")

        return file_groups

    # Populate a list of commands to be run in parallel
    def merge_files(self, in_dir: Path, out_dir: Path) -> list:
        """Function populates a list of
            commands to be executed in parallel
        Args:
            in_dir - Input folder with bclconvert
                        output
            out_dir - Output folder[/FASTQ]
        """
        commands = []
        ext = "001.fastq.gz"  # default ext for fastq files
        file_groups = self.group_fastqs(in_dir, out_dir)
        # pprint(file_groups) # Uncomment this to see fastq file groups
        for merged_file, fastqs in file_groups.items():
            fastq_files = " ".join(fastqs)
            commands.append(f"cat {fastq_files} >> {out_dir}/{merged_file}")

        for file in sorted(os.listdir(in_dir)):
            if not file.endswith(ext):
                commands.append(f"mv {in_dir}/{file} {out_dir}")

        return commands

    # Run the list iteratable in parallel with Popen
    def merge_fastq(self, in_dir: Path, out_dir: Path) -> None:
        """Function runs an iteratable in parallel
        Args:
            commands - [list of commands to run in parallel]
        """
        commands = self.merge_files(in_dir, out_dir)
        print("\nConcatenating Fastq Files in Parallel")
        processes = [Popen(cmd, shell=True) for cmd in commands]
        # Execute the `cmd` (child process) in a new process
        for p in processes:
            # wait for the child process to terminate
            p.wait()

        # Remove the /tmp folder once the merging is complete
        shutil.rmtree(path=in_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge Fastq files.")
    parser.add_argument(
        "-i",
        "--input_path",
        help="Input folder path - (should be inside the in_dir)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        help="Output folder path - (should be inside the in_dir)",
        required=True,
    )
    args = parser.parse_args()

    m = MergeLanes(Path(args.input_path), Path(args.output_path))
