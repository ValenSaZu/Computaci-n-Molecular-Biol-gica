//Probar su implementación en el ejemplo visto en clase: alinear las cadenas AAAC con AGC, verificando que se obtengan todas las soluciones posibles e indicando el score obtenido.

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>

using namespace std;

struct Matrix {
    vector<vector<pair<int, vector<pair<int, int>>>>> Matriz;
    Matrix(int lengthCadena1, int lengthCadena2) {
        Matriz.assign(lengthCadena1 + 1, vector<pair<int, vector<pair<int, int>>>>(lengthCadena2 + 1, { 0, {} }));
    }

    void inicializarMatriz(int lengthCadena1, int lengthCadena2) {
        // Primera fila (j = 0..lengthCadena2)
        for (int j = 1; j <= lengthCadena2; j++) {
            Matriz[0][j].first = (-2) * j;
            Matriz[0][j].second.push_back(make_pair(0, j - 1));
        }

        // Primera columna (i = 0..lengthCadena1)
        for (int i = 1; i <= lengthCadena1; i++) {
            Matriz[i][0].first = (-2) * i;
            Matriz[i][0].second.push_back(make_pair(i - 1, 0));
        }
    }
};

// Llena la matriz
Matrix construccionMatriz(vector<char> cadena1, vector<char> cadena2) {
    Matrix matriz(cadena1.size(), cadena2.size());
    matriz.inicializarMatriz(cadena1.size(), cadena2.size());

    int lengthCadena1 = cadena1.size();
    int lengthCadena2 = cadena2.size();

    const int MATCH = 1;
    const int MISMATCH = -1;
    const int GAP = -2;

    for (int i = 1; i <= lengthCadena1; i++) {
        for (int j = 1; j <= lengthCadena2; j++) {
            int scoreDiagonal = matriz.Matriz[i - 1][j - 1].first +
                (cadena1[i - 1] == cadena2[j - 1] ? MATCH : MISMATCH);
            int scoreArriba = matriz.Matriz[i - 1][j].first + GAP;
            int scoreIzquierda = matriz.Matriz[i][j - 1].first + GAP;

            int mejorScore = max({ scoreDiagonal, scoreArriba, scoreIzquierda });
            matriz.Matriz[i][j].first = mejorScore;

            if (scoreDiagonal == mejorScore) {
                matriz.Matriz[i][j].second.push_back(make_pair(i - 1, j - 1));
            }
            if (scoreArriba == mejorScore) {
                matriz.Matriz[i][j].second.push_back(make_pair(i - 1, j));
            }
            if (scoreIzquierda == mejorScore) {
                matriz.Matriz[i][j].second.push_back(make_pair(i, j - 1));
            }
        }
    }
    return matriz;
}

// Reconstruccion de UN solo alineamiento (sigue siempre el primer predecesor)
pair<string, string> reconstruirAlineamiento(Matrix& matriz, vector<char> cadena1, vector<char> cadena2) {
    string alin1 = "";
    string alin2 = "";

    int i = cadena1.size();
    int j = cadena2.size();

    while (i > 0 || j > 0) {
        pair<int, int> anterior = matriz.Matriz[i][j].second[0];
        int pi = anterior.first;
        int pj = anterior.second;

        if (pi == i - 1 && pj == j - 1) {
            alin1 = cadena1[i - 1] + alin1;
            alin2 = cadena2[j - 1] + alin2;
        }
        else if (pi == i - 1 && pj == j) {
            alin1 = cadena1[i - 1] + alin1;
            alin2 = '-' + alin2;
        }
        else {
            alin1 = '-' + alin1;
            alin2 = cadena2[j - 1] + alin2;
        }

        i = pi;
        j = pj;
    }

    return make_pair(alin1, alin2);
}

// Reconstruccion de todos los caminos
void reconstruirAlineamientos(Matrix& matriz, vector<char>& cadena1, vector<char>& cadena2,
    int i, int j, string alin1, string alin2,
    vector<pair<string, string>>& resultados) {
    if (i == 0 && j == 0) {
        resultados.push_back(make_pair(alin1, alin2));
        return;
    }

    for (auto& anterior : matriz.Matriz[i][j].second) {
        int pi = anterior.first;
        int pj = anterior.second;

        string nuevoAlin1 = alin1;
        string nuevoAlin2 = alin2;

        if (pi == i - 1 && pj == j - 1) {
            nuevoAlin1 = cadena1[i - 1] + nuevoAlin1;
            nuevoAlin2 = cadena2[j - 1] + nuevoAlin2;
        }
        else if (pi == i - 1 && pj == j) {
            nuevoAlin1 = cadena1[i - 1] + nuevoAlin1;
            nuevoAlin2 = '-' + nuevoAlin2;
        }
        else { // pi == i, pj == j - 1
            nuevoAlin1 = '-' + nuevoAlin1;
            nuevoAlin2 = cadena2[j - 1] + nuevoAlin2;
        }

        reconstruirAlineamientos(matriz, cadena1, cadena2, pi, pj, nuevoAlin1, nuevoAlin2, resultados);
    }
}

int main() {
    // Ejemplo: AAAC vs AGC
    vector<char> cadena1 = { 'A', 'A', 'A', 'C' };
    vector<char> cadena2 = { 'A', 'G', 'C' };

    Matrix matriz = construccionMatriz(cadena1, cadena2);

    int scoreOptimo = matriz.Matriz[cadena1.size()][cadena2.size()].first;
    cout << "Cadena 1: AAAC" << endl;
    cout << "Cadena 2: AGC" << endl;
    cout << "Score optimo: " << scoreOptimo << endl << endl;

    vector<pair<string, string>> todosLosAlineamientos;
    reconstruirAlineamientos(matriz, cadena1, cadena2, cadena1.size(), cadena2.size(), "", "", todosLosAlineamientos);

    cout << "Numero de alineamientos optimos encontrados: " << todosLosAlineamientos.size() << endl << endl;

    for (size_t k = 0; k < todosLosAlineamientos.size(); k++) {
        cout << "Alineamiento " << (k + 1) << ":" << endl;
        cout << todosLosAlineamientos[k].first << endl;
        cout << todosLosAlineamientos[k].second << endl;
        cout << endl;
    }

    return 0;
}