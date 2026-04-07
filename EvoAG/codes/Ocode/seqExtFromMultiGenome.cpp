#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

void read_genome_file(const std::string &filename, std::string &genome_seq) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open " << filename << "!\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            continue;
        } else {
            genome_seq += line;
        }
    }
    file.close();
}

void read_positions_file(const std::string &filename, std::string &ag, const std::unordered_map<std::string, std::string> &genomes) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open " << filename << "!\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string sp;
        int st, en;
        if (iss >> sp >> st >> en) {
            int seq_len = en - st + 1;
            st = st - 1;
            auto it = genomes.find(sp);
            if (it != genomes.end()) {
                ag += it->second.substr(st, seq_len);
            } else {
                std::cerr << "Species " << sp << " not found in genome data\n";
            }
        }
    }
    file.close();
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <genome1> <genome2> <positions>\n";
        return EXIT_FAILURE;
    }

    std::unordered_map<std::string, std::string> genomes;
    std::string genome1_file = std::string(argv[1]) + ".fasta";
    std::string genome2_file = std::string(argv[2]) + ".fasta";

    read_genome_file(genome1_file, genomes[argv[1]]);
    read_genome_file(genome2_file, genomes[argv[2]]);

    std::string ag;
    read_positions_file(argv[3], ag, genomes);

    std::cout << ">Salmonella orthologous_ancient genome\n";
    for (size_t i = 0; i < ag.size(); i += 80) {
        std::cout << ag.substr(i, 80) << "\n";
    }

    return EXIT_SUCCESS;
}