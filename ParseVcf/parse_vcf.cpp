#include <fstream>
#include <string>
#include <vector>
#include <iostream>

/*
 *@brief Parse the Variant Call format
 */
class Vcf
{
public:
    std::string chrom;
    std::string pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::string format;
    std::string sample;
    std::vector<std::string> headers;
    std::vector<std::string> variants;
    Vcf();
    Vcf(std::string path);
    ~Vcf();
    std::string getFileFormat();
    std::string getFileDate();
    std::string getFilter();
    std::string getSource();
    std::string getReference();
    std::string getVariants();
    std::string getHeaders();
    void PrintHeaders();
    void PrintVariants();
    void ReadVcf(std::string path);
    void WriteHeaders(std::string path);
    void WriteVariants(std::string path);
    friend std::ofstream &operator<<(std::ofstream &ofs, Vcf &vcf);
    friend std::ifstream &operator>>(std::ifstream &ifs, Vcf &vcf);
    friend std::ostream &operator<<(std::ostream &out, Vcf &vcf);
};
/*
 *@brief Parse a vcf file - Default Constructor
 */
Vcf::Vcf()
{
    std::cout << "Default Constructor Called" << std::endl;
}

/*
 *@brief Parse a vcf file with the given path - Parameterized Constructor
 */
Vcf::Vcf(std::string path)
{
    std::fstream infile(path);
    if (infile.is_open())
    {
        std::string entry;
        while (getline(infile, entry, '\n'))
        {
            if (entry.rfind("##", 0) == 0)
            {
                headers.push_back(entry);
            }
            else
            {
                variants.push_back(entry);
            }
        }
    }
    else
    {
        std::cout << "Failed to open VCF file. Check if it exists." << std::endl;
    }
    for (int i = 0; i < 10; ++i)
    {
        std::cout << variants[i] << std::endl;
    }
    infile.close();
    variants.shrink_to_fit();
    headers.shrink_to_fit();
}

void Vcf::ReadVcf(std::string path)
{
    std::fstream infile(path);
    if (infile.is_open())
    {
        std::string entry;
        while (getline(infile, entry, '\n'))
        {
            if (entry.rfind("##", 0) == 0)
            {
                headers.push_back(entry);
            }
            else
            {
                variants.push_back(entry);
            }
        }
    }
    else
    {
        std::cout << "Failed to open VCF file. Check if it exists." << std::endl;
    }
    infile.close();
    variants.shrink_to_fit();
    headers.shrink_to_fit();
}

/*
 *@brief Get vcf file format version
 */
std::string Vcf::getFileFormat()
{
    std::string fileformat;
    std::string notfound = "File format not found";
    for (std::string line : headers)
    {
        if (line.rfind("##fileformat", 0) == 0)
        {
            int pos = line.find("=");
            fileformat = line.substr(pos + 1);
            return fileformat;
        }
    }
    return notfound;
}

/*
 *@brief Get vcf file date
 */
std::string Vcf::getFileDate()
{
    std::string fileDate;
    std::string notfound = "File Date not found";
    for (std::string line : headers)
    {
        if (line.rfind("##fileDate", 0) == 0)
        {
            int pos = line.find("=");
            fileDate = line.substr(pos + 1);
            return fileDate;
        }
    }
    return notfound;
}

/*
 *@brief Get vcf file filter
 */
std::string Vcf::getFilter()
{
    std::string filter;
    std::string notfound = "Filtering information not found";
    for (std::string line : headers)
    {
        if (line.rfind("##FILTER", 0) == 0)
        {
            int pos = line.find("=");
            filter = line.substr(pos + 1);
            return filter;
        }
    }
    return notfound;
}

/*
 *@brief Get vcf file source
 */
std::string Vcf::getSource()
{
    std::string source;
    std::string notfound = "VCF source not found";
    for (std::string line : headers)
    {
        if (line.rfind("##source", 0) == 0)
        {
            int pos = line.find("=");
            source = line.substr(pos + 1);
            return source;
        }
    }
    return notfound;
}

/*
 *@brief Get vcf file generating reference
 */
std::string Vcf::getReference()
{
    std::string reference;
    std::string notfound = "Reference genome not found";
    for (std::string line : headers)
    {
        if (line.rfind("##reference", 0) == 0)
        {
            int pos = line.find("=");
            reference = line.substr(pos + 1);
            return reference;
        }
    }
    return notfound;
}

/*
 *@brief Get all variants from the vcf file
 */
void Vcf::PrintVariants()
{
    std::string notfound = "Variant(s) not found";
    if (variants.size() > 0)
    {
        for (auto vcf : variants)
        {
            std::cout << vcf << std::endl;
        }
    }
    else
    {
        std::cout << "Variant(s) not found" << std::endl;
    }
}

/*
 *@brief Get all headers from the vcf file
 */
void Vcf::PrintHeaders()
{
    if (headers.size() > 0)
    {
        for (std::string header : headers)
        {
            std::cout << header << std::endl;
        }
    }
    else
    {
        std::cout << "Header(s) not found" << std::endl;
    }
}

/*
 *@brief Overloaded operator to write a vcf object to a file
 */
std::ofstream &operator<<(std::ofstream &ofs, Vcf &vcf)
{
    ofs << vcf.chrom << "\t";
    ofs << vcf.pos << "\t";
    ofs << vcf.id << "\t";
    ofs << vcf.ref << "\t";
    ofs << vcf.alt << "\t";
    ofs << vcf.qual << "\t";
    ofs << vcf.filter << "\t";
    ofs << vcf.info << "\t";
    ofs << vcf.format << "\t";
    ofs << vcf.sample << std::endl;
    return ofs;
}

/*
 *@brief Overloaded operator to retrieve a vcf object from a file
 */
std::ifstream &operator>>(std::ifstream &ifs, Vcf &vcf)
{
    ifs >> vcf.chrom;
    ifs >> vcf.pos;
    ifs >> vcf.id;
    ifs >> vcf.ref;
    ifs >> vcf.alt;
    ifs >> vcf.qual;
    ifs >> vcf.filter;
    ifs >> vcf.info;
    ifs >> vcf.format;
    ifs >> vcf.sample;
    return ifs;
}

/*
 *@brief Overloaded operator to write a vcf object to the console
 */
std::ostream &operator<<(std::ostream &out, Vcf &vcf)
{
    out << vcf.chrom << "\t";
    out << vcf.pos << "\t";
    out << vcf.id << "\t";
    out << vcf.ref << "\t";
    out << vcf.alt << "\t";
    out << vcf.qual << "\t";
    out << vcf.filter << "\t";
    out << vcf.info << "\t";
    out << vcf.format << "\t";
    out << vcf.sample << std::endl;
    return out;
}

/*
*@brief Clear the contents in the headers/variants vectors to
        free up space in heap. - Destructor
*/
Vcf::~Vcf()
{
    headers.clear();
    variants.clear();
}

int main()
{
    Vcf vcf;
    vcf.ReadVcf("example.vcf");
    std::cout << "File Format      : " << vcf.getFileFormat() << std::endl;
    std::cout << "File Date        : " << vcf.getFileDate() << std::endl;
    std::cout << "File Filter      : " << vcf.getFilter() << std::endl;
    std::cout << "VCF Source       : " << vcf.getSource() << std::endl;
    std::cout << "Reference Genome : " << vcf.getReference() << std::endl;
    std::cout << "Variant(s) Found : " << std::endl;
    vcf.PrintVariants();
    return 0;
}
