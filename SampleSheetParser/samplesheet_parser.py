import argparse
import pathlib
from typing import Dict
from pathlib import Path

# from pprint import pprint


class SampleSheet:
    def __init__(self, samplesheet: Path) -> None:
        self.samplesheet = samplesheet
        self.Header: dict = {}
        self.Reads: list = []
        self.Settings: dict = {}
        self.Data: list = []
        # will be updated by getSectionsBoundary()
        self.Sections: dict = {
            "[Header]": {"start": -1, "end": 0},
            "[Settings]": {"start": -1, "end": 0},
            "[Reads]": {"start": -1, "end": 0},
            "[Data]": {"start": -1, "end": 0},
        }
        self.parse(samplesheet)
        # pprint(self.Header)
        # pprint(self.Reads)
        # pprint(self.Settings)
        # pprint(self.Data)

    # Parse the Header Section from the samplesheet
    # Header Section will be returned as a dictionary
    def parse_Header(self, lines: list) -> None:
        """Parse `[Header]` from Illumina SampleSheet
        Args:
            lines - list of all lines in the 
                    samplesheet
        """
        start = self.Sections["[Header]"]["start"]
        end = self.Sections["[Header]"]["end"]
        for line in lines[start:end]:
            if not line:
                continue
            var = line.split(",")
            key, value = var[0], var[1:]
            self.Header[key] = value[0]

    # Parse the Reads Section from the samplesheet
    # Reads Section will be returned as a list
    def parse_Reads(self, lines: list) -> None:
        """Parse `[Reads]` from Illumina SampleSheet
        Args:
            lines - list of all lines in the 
                    samplesheet
        """
        start = self.Sections["[Reads]"]["start"]
        end = self.Sections["[Reads]"]["end"]
        for line in lines[start:end]:
            if not line:
                continue
            key = line.split(",")[0]
            value = line.split(",")[1]
            if value == "":
                value = " "
            self.Reads.append(":".join([key, value]))

    # Parse the Settings Section from the samplesheet
    # Settings Section will be returned as a dictionary
    def parse_Settings(self, lines: list) -> None:
        """Parse `[Settings]` from Illumina SampleSheet
        Args:
            lines - list of all lines in the 
                    samplesheet
        """
        start = self.Sections["[Settings]"]["start"]
        end = self.Sections["[Settings]"]["end"]
        for line in lines[start:end]:
            if not line:
                continue
            var = line.split(",")
            key, value = var[0], var[1:]
            self.Settings[key] = value[0]

    # Parse the Data Section from the samplesheet
    # [Data] Section will be updated in the Data list
    def parse_Data(self, lines: list) -> None:
        """Parse `[Data]` from Illumina SampleSheet
        Args:
            lines - list of all lines in the 
                    samplesheet
        """
        start = self.Sections["[Data]"]["start"]
        end = self.Sections["[Data]"]["end"] + 1
        header = lines[start].strip().split(",")
        for line in lines[1 + start : end]:
            if not line:
                continue
            row: Dict = {}
            for j, value in enumerate(line.strip().split(",")):
                # If the data line is longer than the header
                # terminate the loop
                if j >= len(header):
                    break
                row[header[j]] = value
            self.Data.append(row)

    # Get Boundaries for each Section in the samplesheet
    # getSectionBoundary() updates the sections dictionary
    def getSectionBoundary(self, lines: list) -> None:
        """Get Section Boundaries for all 4 parts of the
        SampleSheet - [Header], [Reads], [Settings],
                      [Data]
        Args:
            lines - list of all lines in the 
                    samplesheet
        """
        i = -1
        curSection = None
        for line in lines:
            i = i + 1
            if not line:
                continue
            var = line.split(",")
            if var[0] in self.Sections.keys():
                curSection = var[0]
                self.Sections[curSection]["start"] = i + 1
            else:
                if curSection:
                    self.Sections[curSection]["end"] = i

    # Parses the Ilumina samplesheet
    # Returns a list of all the lines in the samplesheet
    def parse(self, samplesheet: Path) -> None:
        """Parse Illumina SampleSheet
        Args:
            samplesheet - Path to Illumina Samplesheet
        """
        file = pathlib.Path(samplesheet)
        if file.exists():
            try:
                with open(samplesheet, mode="r") as sheet:
                    lines = sheet.readlines()
                    self.getSectionBoundary(lines)
                    self.parse_Header(lines)
                    self.parse_Reads(lines)
                    self.parse_Settings(lines)
                    self.parse_Data(lines)
            except:
                print("Error Parsing SampleSheet!")
        else:
            print("SampleSheet Path Not Found!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Illumina Samplesheet")
    parser.add_argument(
        "-i", "--samplesheet", help="Input Samplesheet Path", required=True,
    )
    args = parser.parse_args()
    s = SampleSheet((args.samplesheet))
