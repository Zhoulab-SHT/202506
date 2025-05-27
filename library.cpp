#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <vector>

int main() {
    // Path to the input text file
    std::string input_file_path = "input.txt";

    // Path to the output text file
    std::string output_file_path = "output.txt";

    // Open the input file
    std::ifstream input_file(input_file_path);
    if (!input_file.is_open()) {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }

    // Open the output file
    std::ofstream output_file(output_file_path);
    if (!output_file.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return 1;
    }

    // Process each line in the input file
    std::string line;
    while (std::getline(input_file, line)) {
        // Split the line into protein ID and sequence
        std::istringstream iss(line);
        std::string protein_id, sequence;
        if (!(iss >> protein_id >> sequence)) {
            std::cerr << "Invalid line format: " << line << std::endl;
            continue;
        }

        // Extract lowercase fragments separated by uppercase letters
        std::vector<std::string> fragments;
        std::regex uppercase_regex("[A-Z]");
        std::sregex_token_iterator iter(sequence.begin(), sequence.end(), uppercase_regex, -1);
        std::sregex_token_iterator end;

        while (iter != end) {
            std::string fragment = *iter++;
            // Convert fragment to lowercase
            std::transform(fragment.begin(), fragment.end(), fragment.begin(), ::tolower);
            fragments.push_back(fragment);
        }

        // Slice each fragment into 40-character segments
        for (const auto& fragment : fragments) {
            if (fragment.length() > 40) {
                for (size_t i = 0; i < fragment.length(); i += 40) {
                    std::string segment = fragment.substr(i, 40);
                    if (segment.length() == 40) {
                        output_file << protein_id << "\t" << segment << std::endl;
                    }
                }
            }
        }
    }

    // Close the files
    input_file.close();
    output_file.close();

    std::cout << "Processing completed. Results written to " << output_file_path << std::endl;

    return 0;
}