#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <cctype>
#include <sstream>
#include <limits>
#include <chrono>
#include <ctime>

using namespace std;

const int MATCH_SCORE = 1;
const int MISMATCH_PENALTY = -1;
const int GAP_PENALTY = -2;


// ============================================================
// DATA STRUCTURES
// ============================================================

struct Sequence
{
    string id;
    string species;
    string type;
    string sequence;
};


struct AlignmentResult
{
    string alignedA;
    string alignedB;
    string symbols;

    int score = 0;

    int matches = 0;
    int mismatches = 0;
    int gaps = 0;

    double identity = 0.0;
};


struct DistanceMatrix
{
    vector<string> labels;
    vector<vector<double>> values;
};

struct PhyloNode
{
    string label;
    PhyloNode* left = nullptr;
    PhyloNode* right = nullptr;
    double leftLength = 0.0;
    double rightLength = 0.0;
};


// ============================================================
// FASTA READER
// ============================================================

Sequence readFASTA(
    const string& filename,
    const string& species,
    const string& type)
{
    ifstream file(filename);

    Sequence result;

    result.species = species;
    result.type = type;

    if (!file.is_open())
    {
        cerr << "ERROR: Could not open "
             << filename << endl;

        return result;
    }

    string line;

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        if (line[0] == '>')
        {
            result.id = line.substr(1);
            continue;
        }

        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        result.sequence += line;
    }

    file.close();

    transform(
        result.sequence.begin(),
        result.sequence.end(),
        result.sequence.begin(),
        [](unsigned char c)
        {
            return static_cast<char>(toupper(c));
        }
    );

    return result;
}


// ============================================================
// SEQUENCE VALIDATION
// ============================================================

bool validateSequence(const Sequence& seq)
{
    if (seq.sequence.empty())
    {
        cerr << "Eroare: Secventa fara nucleotide: "
             << seq.species << endl;

        return false;
    }

    for (char c : seq.sequence)
    {
        if (c != 'A' &&
            c != 'C' &&
            c != 'G' &&
            c != 'T' &&
            c != 'N')
        {
            cerr << "Eroare:Secventa invalida '"
                 << c
                 << "' in "
                 << seq.species
                 << endl;

            return false;
        }
    }

    return true;
}


// ============================================================
// NEEDLEMAN-WUNSCH
// ============================================================

AlignmentResult needlemanWunsch(
    const string& seqA,
    const string& seqB)
{
    int n = static_cast<int>(seqA.length());
    int m = static_cast<int>(seqB.length());

    vector<vector<int>> F(
        n + 1,
        vector<int>(m + 1, 0)
    );


    // --------------------------------------------------------
    // INITIALIZATION
    // --------------------------------------------------------

    for (int i = 0; i <= n; ++i)
        F[i][0] = i * GAP_PENALTY;

    for (int j = 0; j <= m; ++j)
        F[0][j] = j * GAP_PENALTY;


    // --------------------------------------------------------
    // DYNAMIC PROGRAMMING
    // --------------------------------------------------------

    for (int i = 1; i <= n; ++i)
    {
        for (int j = 1; j <= m; ++j)
        {
            int substitution =
                (seqA[i - 1] == seqB[j - 1])
                ? MATCH_SCORE
                : MISMATCH_PENALTY;

            int diagonal =
                F[i - 1][j - 1] + substitution;

            int up =
                F[i - 1][j] + GAP_PENALTY;

            int left =
                F[i][j - 1] + GAP_PENALTY;

            F[i][j] =
                max({
                    diagonal,
                    up,
                    left
                });
        }
    }


    // --------------------------------------------------------
    // TRACEBACK
    // --------------------------------------------------------

    string alignA;
    string alignB;
    string symbols;

    int i = n;
    int j = m;

    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0)
        {
            int substitution =
                (seqA[i - 1] == seqB[j - 1])
                ? MATCH_SCORE
                : MISMATCH_PENALTY;

            if (F[i][j] ==
                F[i - 1][j - 1] + substitution)
            {
                alignA += seqA[i - 1];
                alignB += seqB[j - 1];

                if (seqA[i - 1] == seqB[j - 1])
                    symbols += '|';
                else
                    symbols += '.';

                --i;
                --j;

                continue;
            }
        }


        if (i > 0 &&
            F[i][j] ==
            F[i - 1][j] + GAP_PENALTY)
        {
            alignA += seqA[i - 1];
            alignB += '-';
            symbols += ' ';

            --i;

            continue;
        }


        alignA += '-';
        alignB += seqB[j - 1];
        symbols += ' ';

        --j;
    }


    reverse(
        alignA.begin(),
        alignA.end()
    );

    reverse(
        alignB.begin(),
        alignB.end()
    );

    reverse(
        symbols.begin(),
        symbols.end()
    );


    // --------------------------------------------------------
    // STATISTICS
    // --------------------------------------------------------

    AlignmentResult result;

    result.alignedA = alignA;
    result.alignedB = alignB;
    result.symbols = symbols;
    result.score = F[n][m];


    for (size_t k = 0;
         k < alignA.length();
         ++k)
    {
        if (alignA[k] == '-' ||
            alignB[k] == '-')
        {
            result.gaps++;
        }
        else if (alignA[k] == alignB[k])
        {
            result.matches++;
        }
        else
        {
            result.mismatches++;
        }
    }


    if (!alignA.empty())
    {
        result.identity =
            100.0 *
            result.matches /
            static_cast<double>(
                alignA.length()
            );
    }

    return result;
}


// ============================================================
// DISTANCE MATRIX
// ============================================================

// The matching fraction is the same quantity currently reported
// as Identity by Needleman-Wunsch: matches / alignment length.
double calculateMatchingFraction(const AlignmentResult& result)
{
    if (result.alignedA.empty())
        return 0.0;

    return static_cast<double>(result.matches) /
           static_cast<double>(result.alignedA.length());
}

// Convert matching fraction into an evolutionary distance.
// d = 1 - matching fraction.
double calculateDistance(const AlignmentResult& result)
{
    return 1.0 - calculateMatchingFraction(result);
}

string makeSequenceLabel(const Sequence& seq)
{
    // Species alone is not always unique (e.g. Homo sapiens appears
    // with both myoglobin and hemoglobin), so the tree uses a label
    // that identifies the actual sequence unambiguously.
    return seq.species + " | " + seq.type;
}

DistanceMatrix buildDistanceMatrix(const vector<Sequence>& sequences)
{
    DistanceMatrix matrix;
    const size_t n = sequences.size();

    matrix.labels.reserve(n);
    for (const auto& seq : sequences)
        matrix.labels.push_back(makeSequenceLabel(seq));

    matrix.values.assign(n, vector<double>(n, 0.0));

    cout << "\n========================================\n";
    cout << "CONSTRUIRE DISTANCE MATRIX\n";
    cout << "========================================\n";
    cout << "Nr perechi:  " << (n * (n - 1)) / 2 << "\n\n";

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            AlignmentResult result = needlemanWunsch(
                sequences[i].sequence,
                sequences[j].sequence
            );

            double distance = calculateDistance(result);
            matrix.values[i][j] = distance;
            matrix.values[j][i] = distance;
        }
    }

    return matrix;
}

void printDistanceMatrix(const DistanceMatrix& matrix)
{
    cout << "\n========================================\n";
    cout << "PAIRWISE DISTANCE MATRIX\n";
    cout << "========================================\n\n";

    const int width = 20;

    cout << left << setw(width) << "Sequence";
    for (const auto& label : matrix.labels)
        cout << right << setw(12) << label.substr(0, 10);
    cout << "\n";

    for (size_t i = 0; i < matrix.labels.size(); ++i)
    {
        cout << left << setw(width)
             << matrix.labels[i].substr(0, width - 2);

        for (size_t j = 0; j < matrix.labels.size(); ++j)
        {
            cout << right << setw(12)
                 << fixed << setprecision(4)
                 << matrix.values[i][j];
        }
        cout << "\n";
    }
}


// ============================================================
// NEIGHBOUR JOINING
// ============================================================

PhyloNode* createLeaf(const string& label)
{
    PhyloNode* node = new PhyloNode;
    node->label = label;
    return node;
}

PhyloNode* createInternalNode(
    PhyloNode* left,
    PhyloNode* right,
    double leftLength,
    double rightLength)
{
    PhyloNode* node = new PhyloNode;
    node->left = left;
    node->right = right;
    node->leftLength = max(0.0, leftLength);
    node->rightLength = max(0.0, rightLength);
    return node;
}

void deletePhyloTree(PhyloNode* node)
{
    if (!node)
        return;

    deletePhyloTree(node->left);
    deletePhyloTree(node->right);
    delete node;
}

string formatBranchLength(double length)
{
    ostringstream out;
    out << fixed << setprecision(4) << max(0.0, length);
    return out.str();
}

string toNewick(const PhyloNode* node)
{
    if (!node)
        return "";

    if (!node->left && !node->right)
        return node->label;

    string left = toNewick(node->left) + ":" +
                  formatBranchLength(node->leftLength);
    string right = toNewick(node->right) + ":" +
                   formatBranchLength(node->rightLength);

    return "(" + left + "," + right + ")";
}

string makeNewickSafeLabel(string label)
{
    // Keep the biological identity readable while avoiding punctuation
    // that has structural meaning in Newick format.
    for (char& c : label)
    {
        if (c == ' ' || c == '|' || c == ':' || c == ',' ||
            c == '(' || c == ')' || c == ';')
        {
            c = '_';
        }
    }
    return label;
}

string toNewickSafe(const PhyloNode* node)
{
    if (!node)
        return "";

    if (!node->left && !node->right)
        return makeNewickSafeLabel(node->label);

    string left = toNewickSafe(node->left) + ":" +
                  formatBranchLength(node->leftLength);
    string right = toNewickSafe(node->right) + ":" +
                   formatBranchLength(node->rightLength);

    return "(" + left + "," + right + ")";
}

string makeTimestamp()
{
    const auto now = chrono::system_clock::now();
    const time_t timeNow = chrono::system_clock::to_time_t(now);
    const auto milliseconds =
        chrono::duration_cast<chrono::milliseconds>(
            now.time_since_epoch()
        ) % 1000;

    tm localTime{};

#ifdef _WIN32
    localtime_s(&localTime, &timeNow);
#else
    localtime_r(&timeNow, &localTime);
#endif

    ostringstream out;
    out << put_time(&localTime, "%Y%m%d_%H%M%S")
        << "_" << setfill('0') << setw(3) << milliseconds.count();
    return out.str();
}

bool saveNewickWithHistory(const PhyloNode* root)
{
    if (!root)
        return false;

    const filesystem::path resultsDir = "results";
    const filesystem::path historyDir = resultsDir / "history";

    error_code ec;
    filesystem::create_directories(historyDir, ec);
    if (ec)
    {
        cerr << "Eroare: Nu s-au putut genera rezultate: "
             << ec.message() << endl;
        return false;
    }

    const string newick = toNewickSafe(root) + ";\n";
    const string timestamp = makeTimestamp();

    const filesystem::path currentFile = resultsDir / "tree.newick";
    const filesystem::path historyFile =
        historyDir / ("tree_" + timestamp + ".newick");

    ofstream current(currentFile);
    ofstream history(historyFile);

    if (!current.is_open() || !history.is_open())
    {
        cerr << "Eroare: Nu s-au putut afisa scorurile Newick" << endl;
        return false;
    }

    current << newick;
    history << newick;

    cout << "\nNewick saved to: " << currentFile.string() << endl;
    cout << "History copy:    " << historyFile.string() << endl;
    return true;
}

PhyloNode* neighborJoining(const DistanceMatrix& input)
{
    const size_t sequenceCount = input.labels.size();

    if (sequenceCount < 2)
    {
        cerr << "ERROR: Neighbor Joining requires at least 2 sequences.\n";
        return nullptr;
    }

    vector<PhyloNode*> nodes;
    vector<vector<double>> distances = input.values;

    for (const auto& label : input.labels)
        nodes.push_back(createLeaf(label));

    while (nodes.size() > 2)
    {
        const size_t n = nodes.size();
        vector<double> rowSums(n, 0.0);

        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                rowSums[i] += distances[i][j];

        double bestQ = numeric_limits<double>::infinity();
        size_t bestI = 0;
        size_t bestJ = 1;

        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = i + 1; j < n; ++j)
            {
                double q =
                    static_cast<double>(n - 2) * distances[i][j]
                    - rowSums[i]
                    - rowSums[j];

                if (q < bestQ)
                {
                    bestQ = q;
                    bestI = i;
                    bestJ = j;
                }
            }
        }

        double dij = distances[bestI][bestJ];

        double leftLength =
            0.5 * dij +
            (rowSums[bestI] - rowSums[bestJ]) /
            (2.0 * static_cast<double>(n - 2));

        double rightLength = dij - leftLength;

        PhyloNode* merged = createInternalNode(
            nodes[bestI], nodes[bestJ], leftLength, rightLength
        );

        vector<PhyloNode*> newNodes;
        vector<size_t> oldIndices;

        for (size_t k = 0; k < n; ++k)
        {
            if (k != bestI && k != bestJ)
            {
                newNodes.push_back(nodes[k]);
                oldIndices.push_back(k);
            }
        }
        newNodes.push_back(merged);

        const size_t newN = newNodes.size();
        vector<vector<double>> newDistances(
            newN, vector<double>(newN, 0.0)
        );

        for (size_t a = 0; a < oldIndices.size(); ++a)
        {
            for (size_t b = 0; b < oldIndices.size(); ++b)
            {
                newDistances[a][b] =
                    distances[oldIndices[a]][oldIndices[b]];
            }
        }

        for (size_t a = 0; a < oldIndices.size(); ++a)
        {
            size_t k = oldIndices[a];
            double newDistance =
                0.5 * (
                    distances[bestI][k] +
                    distances[bestJ][k] -
                    distances[bestI][bestJ]
                );

            newDistances[a][newN - 1] = newDistance;
            newDistances[newN - 1][a] = newDistance;
        }

        nodes.swap(newNodes);
        distances.swap(newDistances);
    }

    double finalDistance = distances[0][1];

    return createInternalNode(
        nodes[0], nodes[1], finalDistance / 2.0, finalDistance / 2.0
    );
}

void printPhylogeneticTree(const PhyloNode* root)
{
    if (!root)
        return;

    cout << "\n========================================\n";
    cout << "ARBORE FILOGENETIC\n";
    cout << "========================================\n\n";
    cout << "Reprezentare Newick:\n";
    cout << toNewick(root) << ";\n";
}

void runPhylogeneticAnalysis(const vector<Sequence>& sequences)
{
    if (sequences.size() < 2)
    {
        cerr << "EROARE: Cel putin 2 secvente valide sunt necesare.\n";
        return;
    }

    // The main sequence vector is the single source of truth.
    // Any sequence added there is automatically included in both
    // the distance matrix and the NJ tree on the next run.
    DistanceMatrix matrix = buildDistanceMatrix(sequences);
    printDistanceMatrix(matrix);

    PhyloNode* root = neighborJoining(matrix);
    printPhylogeneticTree(root);
    saveNewickWithHistory(root);

    deletePhyloTree(root);
}


// ============================================================
// PRINT ALIGNMENT
// ============================================================

void printAlignment(
    const Sequence& A,
    const Sequence& B,
    const AlignmentResult& result)
{
    cout << "\n========================================\n";

    cout << A.species
         << " vs "
         << B.species
         << "\n";

    cout << "========================================\n";

    cout << "Type A:     "
         << A.type
         << "\n";

    cout << "Type B:     "
         << B.type
         << "\n\n";

    cout << "Scor:      "
         << result.score
         << "\n";

    cout << "Similaritate:   "
         << fixed
         << setprecision(2)
         << result.identity
         << "%\n";

    cout << "Matches:    "
         << result.matches
         << "\n";

    cout << "Mismatches: "
         << result.mismatches
         << "\n";

    cout << "Gaps:       "
         << result.gaps
         << "\n\n";


    const int blockSize = 60;

    for (size_t k = 0;
         k < result.alignedA.length();
         k += blockSize)
    {
        cout << A.species
             << ": "
             << result.alignedA.substr(
                    k,
                    blockSize)
             << "\n";

        cout << "          "
             << result.symbols.substr(
                    k,
                    blockSize)
             << "\n";

        cout << B.species
             << ": "
             << result.alignedB.substr(
                    k,
                    blockSize)
             << "\n\n";
    }
}


// ============================================================
// RUN ONE PAIR
// ============================================================

void runPair(
    const vector<Sequence>& sequences,
    int a,
    int b)
{
    const Sequence& seqA =
        sequences[a];

    const Sequence& seqB =
        sequences[b];

    cout << "\nRulare Needleman-Wunsch...\n";

    AlignmentResult result =
        needlemanWunsch(
            seqA.sequence,
            seqB.sequence
        );

    printAlignment(
        seqA,
        seqB,
        result
    );
}


// ============================================================
// RUN ALL PAIRS
// ============================================================

void runAllPairs(
    const vector<Sequence>& sequences)
{
    int totalPairs = 0;

    for (size_t i = 0;
         i < sequences.size();
         ++i)
    {
        for (size_t j = i + 1;
             j < sequences.size();
             ++j)
        {
            totalPairs++;
        }
    }

    cout << "\n========================================\n";
    cout << "ALINIEREA TUTUROR PERECHILOR POSIBILE\n";
    cout << "========================================\n";

    cout << "Nr de secvente: "
         << sequences.size()
         << "\n";

    cout << "Nr de perechi: "
         << totalPairs
         << "\n\n";


    int pairNumber = 0;

    for (size_t i = 0;
         i < sequences.size();
         ++i)
    {
        for (size_t j = i + 1;
             j < sequences.size();
             ++j)
        {
            pairNumber++;

            cout << "\n----------------------------------------\n";

            cout << "Pair "
                 << pairNumber
                 << " / "
                 << totalPairs
                 << "\n";

            cout << sequences[i].species
                 << " vs "
                 << sequences[j].species
                 << "\n";

            cout << "----------------------------------------\n";


            AlignmentResult result =
                needlemanWunsch(
                    sequences[i].sequence,
                    sequences[j].sequence
                );


            cout << "Scor:    "
                 << result.score
                 << "\n";

            cout << "Similaritate: "
                 << fixed
                 << setprecision(2)
                 << result.identity
                 << "%\n";

            cout << "Matches:  "
                 << result.matches
                 << "\n";

            cout << "Mismatches: "
                 << result.mismatches
                 << "\n";

            cout << "Gaps:     "
                 << result.gaps
                 << "\n";
        }
    }

    cout << "\n========================================\n";
    cout << "TOATE PERECHILE POSIBILE AU FOST ALINIATE\n";
    cout << "========================================\n";
}


// ============================================================
// MAIN
// ============================================================

int main()
{
    vector<Sequence> sequences;


    // --------------------------------------------------------
    // SEQUENCE DATABASE
    // --------------------------------------------------------

    sequences.push_back(
        readFASTA(
            "data/mycobacterium.fna",
            "Mycobacterium",
            "Globin-associated sensor"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/bacillus_flavohemoglobin.fna",
            "Bacillus subtilis",
            "Flavohemoglobin"
        )
    );
    sequences.push_back(
        readFASTA(
            "data/saccharomyces_flavohemoglobin.fna",
            "Saccharomyces cerivisiae",
            "Flavohemoglobin"
        )
    );
    sequences.push_back(
        readFASTA(
            "data/lotus_phytoglobin.fna",
            "Lotus japonicus",
            "Nonsymbiotic phytoglobin"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/arabidopsis_phytoglobin.fna",
            "Arabidopsis thaliana",
            "Nonsymbiotic phytoglobin 1"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/caenorhabditis_globin.fna",
            "Caenorhabditis elegans",
            "Globin-like 9"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/mus_cytoglobin.fna",
            "Mus musculus",
            "Cytoglobin"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/anopheles_myoglobin.fna",
            "Anopheles",
            "Myoglobin"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/human_myoglobin.fna",
            "Homo sapiens",
            "Myoglobin"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/xenopus_hemoglobin.fna",
            "Xenopus",
            "Hemoglobin"
        )
    );

    sequences.push_back(
        readFASTA(
            "data/human_hemoglobin.fna",
            "Homo sapiens",
            "Hemoglobin"
        )
    );


    // --------------------------------------------------------
    // VALIDATION
    // --------------------------------------------------------

    cout << "\n========================================\n";
    cout << "       SECVENTELE NUCLEOTIDICE CE CODIFICA PROTEINELE DIN F. GLOBINELOR \n";
    cout << "========================================\n\n";


    int validSequences = 0;

    for (const auto& seq : sequences)
    {
        bool valid =
            validateSequence(seq);

        cout << seq.species
             << " | "
             << seq.type
             << " | "
             << seq.sequence.length()
             << " nt";

        if (valid)
        {
            cout << " | VALID";
            validSequences++;
        }
        else
        {
            cout << " | INVALID";
        }

        cout << "\n";
    }


    cout << "\nSecventele sunt valide: "
         << validSequences
         << " / "
         << sequences.size()
         << "\n";


    if (validSequences !=
        static_cast<int>(sequences.size()))
    {
        cerr << "\nEroare: Unele secvente sunt invalide.\n";
        return 1;
    }


    // --------------------------------------------------------
    // INTERACTIVE MENU
    // --------------------------------------------------------

    while (true)
    {
        cout << "\n\n========================================\n";
        cout << "   ANALIZA SECVENTELOR NUCLEOTIDICE CE CODIFICA PROTEINE DIN F. GLOBINELOR\n";
        cout << "========================================\n\n";

        cout << "Secvente nucleotidice disponibile:\n\n";


        for (size_t i = 0;
             i < sequences.size();
             ++i)
        {
            cout << "["
                 << i + 1
                 << "] "
                 << sequences[i].species
                 << " - "
                 << sequences[i].type
                 << " ("
                 << sequences[i].sequence.length()
                 << " nt)"
                 << "\n";
        }


        cout << "\nOptiuni:\n";
        cout << "  P = Alinierea unei perechi la alegere\n";
        cout << "  T = Alinierea tuturor perechilor posibile\n";
        cout << "  A = Construire arbore filogenetic (NJ)\n";
        cout << "  I = Incheie\n";


        char option;

        cout << "\nSelecteaza optiune: ";
        cin >> option;

        option =
            static_cast<char>(
                toupper(
                    static_cast<unsigned char>(option)
                )
            );


        // ----------------------------------------------------
        // QUIT
        // ----------------------------------------------------

        if (option == 'I')
        {
            cout << "\nProgram incheiat.\n";
            break;
        }


        // ----------------------------------------------------
        // ALL PAIRS
        // ----------------------------------------------------

        if (option == 'T')
        {
            runAllPairs(sequences);
            continue;
        }


        // ----------------------------------------------------
        // PHYLOGENETIC TREE
        // ----------------------------------------------------

        if (option == 'A')
        {
            runPhylogeneticAnalysis(sequences);
            continue;
        }


        // ----------------------------------------------------
        // PAIRWISE
        // ----------------------------------------------------

        if (option == 'P')
        {
            int a;
            int b;


            cout << "\nSelecteaza prima secventa: ";
            cin >> a;

            cout << "Selecteaza a doua secventa: ";
            cin >> b;


            if (a < 1 ||
                b < 1 ||
                a > static_cast<int>(sequences.size()) ||
                b > static_cast<int>(sequences.size()) ||
                a == b)
            {
                cerr << "\nSelectie invalida.\n";
                continue;
            }


            runPair(
                sequences,
                a - 1,
                b - 1
            );

            continue;
        }


        cout << "\nAceasta nu este o optiune.\n";
    }


    return 0;
}