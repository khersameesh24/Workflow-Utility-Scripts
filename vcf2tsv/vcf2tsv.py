import io
import argparse
import pandas as pd


# Class which contains utility functions to convert VCF to TSV
class VcfToTsv:
    """
    Read a vcf file and convert it to a tsv file.

    Uses:
        `Pandas`
    """

    def __init__(self, input_vcf: str, output_tsv: str):
        self.input_vcf = input_vcf
        self.output_tsv = output_tsv

        # Read VCF
        lines: list = self.read_vcf()

        # Write output to TSV file
        self.write_output(lines)

    # Read data in the input VCF file
    def read_vcf(self) -> list:
        try:
            with open(self.input_vcf, "rt") as f:
                lines: list = [l for l in f if not l.startswith("##")]
                return lines
        except Exception:
            print("Error Reading VCF File")

    def write_output(self, lines: list) -> None:
        data = pd.read_csv(
            io.StringIO("".join(lines)),
            dtype={
                "#CHROM": str,
                "POS": int,
                "ID": str,
                "REF": str,
                "ALT": str,
                "QUAL": str,
                "FILTER": str,
                "INFO": str,
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})
        main_vcf = data.filter(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"], axis=1
        )
        info_data_modified = self.get_info_data(data)
        format_data_modified = self.get_format_data(data)

        modified_tags = main_vcf.join(info_data_modified)
        modified_tags = modified_tags.join(format_data_modified)
        modified_tags.to_csv(self.output_tsv, sep="\t", encoding="utf-8", index=False)

    # Parse Info Data
    def get_info_data(self, data: pd) -> pd:
        # Get original info data
        info = data["INFO"].str.split(";", expand=True)

        # Get unique tags in Info field and create a dataframe using these tags as columns
        tags: list = []
        for i in info.columns:
            tags.extend(info[i].str.split("=", expand=True)[0].tolist())
        not_none_tags = set(filter(None.__ne__, tags))
        info_data = pd.DataFrame(columns=not_none_tags)

        # Parse info data from the VCF and append data to main dataframe - info_data
        for index, row in info.iterrows():
            # get data for each INFO row and convert that to a dataframe
            entry = row.str.split("=", expand=True).transpose()
            new_header = entry.iloc[0]
            entry = entry[1:]
            entry.columns = new_header
            # Drop empty columns
            entry = entry.loc[:, entry.columns.notnull()].fillna("Yes")
            # Concatenate with the main info dataframe
            info_data = pd.concat([info_data, entry])
        # Reset index and replace NaN with empty string
        info_data = info_data.reset_index(drop=True)
        info_data = info_data.fillna("")

        return info_data

    # Parse FORMAT Data
    def get_format_data(self, data: pd) -> pd:
        # Get the VCF #CHROM line
        header = list(data.head(0))

        # Determine the sample columns which appear after the FORMAT column in the VCF file
        samples = self.get_sample_names(header)

        # Get the tags in format column
        format_tags = data["FORMAT"].str.split(":", expand=True)

        # Construct an overall dataframe from the unique tags in format column
        unique_column_names = pd.unique(format_tags.stack())
        format_tags_unique = pd.DataFrame(columns=unique_column_names)

        for sample in samples:
            format_tags_sample = pd.DataFrame(columns=unique_column_names)
            # Append the sample name to FORMAT tags which is the header in the tsv file
            format_names = [sub + "-" + sample for sub in unique_column_names]
            # Fetch the format data for the sample
            sample_data = data[sample].str.split(":", expand=True)

            for i in sample_data.index:
                # Construct a dataframe for each row in the sample data
                head = format_tags.loc[[i]].iloc[0]
                value = sample_data.loc[[i]]
                value.columns = head

                # Skip columns which have 'None' and reset index
                value = value.loc[:, value.columns.notnull()]
                value = value.reset_index(drop=True)

                # Concatenate the sample data to overall sample dataframe
                format_tags_sample = pd.concat([format_tags_sample, value]).fillna("")
            # Concatenate the complete sample data to the
            # main dataframe containing unique tags and reset index
            format_tags_sample.columns = format_names
            format_tags_sample = format_tags_sample.reset_index(drop=True)
            if format_tags_unique.empty:
                format_tags_unique = format_tags_sample
            else:
                format_tags_unique = format_tags_unique.join(format_tags_sample)

        return format_tags_unique

    # Get Sample Names from VCF file
    def get_sample_names(self, header: list) -> list:
        samples: list = []
        for i in range((header.index("FORMAT") + 1), (len(header))):
            samples.append(header[i])
        return samples


# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert vcf to tsv file")
    parser.add_argument(
        "-i", "--input_vcf", help="Input VCF file", required=True,
    )
    parser.add_argument(
        "-o", "--output_tsv", help="Output TSV file", required=True,
    )
    args = parser.parse_args()
    VcfToTsv(input_vcf=args.input_vcf, output_tsv=args.output_tsv)
