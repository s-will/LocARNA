#include <iostream>
#include <vector>
#include <string>

using namespace std;

typedef string sequence_t;

/************************************************************
 *
 * compute identity of two sequences, i.e.
 * the best identity of any alignment of the two sequences
 *
 * identity of an alignment = number of matchs / length of alignment
 *
 * for computing this identity it suffices to maximize the number of matchs,
 * since in the same time the length of the alignment is minimized
 *
 * USAGE: identity [<seq1> <seq2>]
 * without arguments read from stdin in simplified mfasta format
 *
 ***********************************************************/

double
identity(sequence_t a, sequence_t b) {
    const size_t n = a.length();
    const size_t m = b.length();

    vector<vector<uint> > M;
    M.resize(n + 1);
    for (uint i = 0; i <= n; i++)
        M[i].resize(m + 1);

    for (uint i = 0; i <= n; i++)
        M[i][0] = 0;
    for (uint j = 0; j <= m; j++)
        M[0][j] = 0;

    for (uint i = 1; i <= n; i++)
        for (uint j = 1; j <= m; j++)
            M[i][j] = max(M[i - 1][j - 1] + ((a[i - 1] == b[j - 1]) ? 1 : 0),
                          max(M[i - 1][j], M[i][j - 1]));

    return M[n][m] / (double)(n + m - M[n][m]);
}

int
main(int argc, char **argv) {
    string seq1;
    string seq2;
    if (argc == 3) {
        seq1 = argv[1];
        seq2 = argv[2];
    } else if (argc == 1) {
        string line;
        getline(cin, line);
        if (line[0] == '>')
            getline(cin, line);
        seq1 = line;

        getline(cin, line);
        if (line[0] == '>')
            getline(cin, line);
        seq2 = line;
    } else {
        cerr << "USAGE: identity [<seq1> <seq2>]" << endl;
        exit(-1);
    }

    cout << identity(seq1, seq2) << endl;
    exit(0);
}
