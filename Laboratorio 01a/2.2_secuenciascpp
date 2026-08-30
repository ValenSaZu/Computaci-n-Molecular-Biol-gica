#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <string>
#include <map>
#include <chrono>
#include <cctype>
#include <climits>

using namespace std;

struct Matrix {
    vector<vector<pair<int, vector<pair<int, int>>>>> Matriz;
    Matrix(int lengthCadena1, int lengthCadena2) {
        Matriz.assign(lengthCadena1 + 1, vector<pair<int, vector<pair<int, int>>>>(lengthCadena2 + 1, { 0, {} }));
    }

    void inicializarMatriz(int lengthCadena1, int lengthCadena2) {
        for (int j = 1; j <= lengthCadena2; j++) {
            Matriz[0][j].first = (-2) * j;
            Matriz[0][j].second.push_back(make_pair(0, j - 1));
        }

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

// Reconstruccion de un solo alineamiento
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

bool reconstruirAlineamientos(Matrix& matriz, vector<char>& cadena1, vector<char>& cadena2,
    int i, int j, string alin1, string alin2,
    vector<pair<string, string>>& resultados,
    size_t limiteResultados) {
    if (resultados.size() >= limiteResultados) {
        return true;
    }

    if (i == 0 && j == 0) {
        resultados.push_back(make_pair(alin1, alin2));
        return resultados.size() >= limiteResultados;
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
        else {
            nuevoAlin1 = '-' + nuevoAlin1;
            nuevoAlin2 = cadena2[j - 1] + nuevoAlin2;
        }

        bool limiteAlcanzado = reconstruirAlineamientos(matriz, cadena1, cadena2, pi, pj,
            nuevoAlin1, nuevoAlin2,
            resultados, limiteResultados);
        if (limiteAlcanzado) return true;
    }
    return false;
}

// Cuenta cuantos caminos optimos hay
long long contarCaminos(Matrix& matriz, int i, int j, vector<vector<long long>>& memo) {
    if (i == 0 && j == 0) return 1;
    if (memo[i][j] != -1) return memo[i][j];

    long long total = 0;
    for (auto& anterior : matriz.Matriz[i][j].second) {
        long long sub = contarCaminos(matriz, anterior.first, anterior.second, memo);
        if (total > 0 && sub > (numeric_limits<long long>::max() - total)) {
            total = numeric_limits<long long>::max();
        }
        else {
            total += sub;
        }
    }
    memo[i][j] = total;
    return total;
}

enum TipoMovimiento { NINGUNO = 0, DIAGONAL = 1, ARRIBA = 2, IZQUIERDA = 3 };

struct ResultadoAlineamiento {
    string alin1;
    string alin2;
    int rupturas;
};

ResultadoAlineamiento seleccionarMenosRupturas(Matrix& matriz, vector<char>& cadena1, vector<char>& cadena2) {
    int n = (int)cadena1.size();
    int m = (int)cadena2.size();
    const int INF = INT_MAX / 2;

    vector<vector<array<int, 4>>> costo(n + 1, vector<array<int, 4>>(m + 1, { INF, INF, INF, INF }));
    vector<vector<array<pair<int, int>, 4>>> predCelda(n + 1,
        vector<array<pair<int, int>, 4>>(m + 1, { make_pair(-1,-1), make_pair(-1,-1), make_pair(-1,-1), make_pair(-1,-1) }));
    vector<vector<array<int, 4>>> predTipo(n + 1, vector<array<int, 4>>(m + 1, { NINGUNO, NINGUNO, NINGUNO, NINGUNO }));

    costo[0][0][NINGUNO] = 0;

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            if (i == 0 && j == 0) continue;

            for (auto& anterior : matriz.Matriz[i][j].second) {
                int pi = anterior.first;
                int pj = anterior.second;

                int movActual;
                if (pi == i - 1 && pj == j - 1) movActual = DIAGONAL;
                else if (pi == i - 1 && pj == j) movActual = ARRIBA;
                else movActual = IZQUIERDA;

                for (int tipoPred = 0; tipoPred < 4; tipoPred++) {
                    int costoPred = costo[pi][pj][tipoPred];
                    if (costoPred >= INF) continue;

                    // Se suma una ruptura solo si el movimiento actual es un gap
                    int incremento = (movActual != DIAGONAL && movActual != tipoPred) ? 1 : 0;
                    int candidato = costoPred + incremento;

                    if (candidato < costo[i][j][movActual]) {
                        costo[i][j][movActual] = candidato;
                        predCelda[i][j][movActual] = make_pair(pi, pj);
                        predTipo[i][j][movActual] = tipoPred;
                    }
                }
            }
        }
    }

    int mejorTipo = DIAGONAL;
    for (int t = DIAGONAL; t <= IZQUIERDA; t++) {
        if (costo[n][m][t] < costo[n][m][mejorTipo]) mejorTipo = t;
    }

    ResultadoAlineamiento resultado;
    resultado.rupturas = costo[n][m][mejorTipo];

    string alin1 = "", alin2 = "";
    int i = n, j = m, tipo = mejorTipo;
    while (!(i == 0 && j == 0)) {
        if (tipo == DIAGONAL) {
            alin1 = cadena1[i - 1] + alin1;
            alin2 = cadena2[j - 1] + alin2;
        }
        else if (tipo == ARRIBA) {
            alin1 = cadena1[i - 1] + alin1;
            alin2 = '-' + alin2;
        }
        else { // IZQUIERDA
            alin1 = '-' + alin1;
            alin2 = cadena2[j - 1] + alin2;
        }

        pair<int, int> pred = predCelda[i][j][tipo];
        int tipoPred = predTipo[i][j][tipo];
        i = pred.first;
        j = pred.second;
        tipo = tipoPred;
    }

    resultado.alin1 = alin1;
    resultado.alin2 = alin2;
    return resultado;
}

int contarRupturas(const string& alin1, const string& alin2) {
    int rupturas = 0;
    int tipoAnterior = NINGUNO;

    for (size_t k = 0; k < alin1.size(); k++) {
        int tipoActual;
        if (alin1[k] == '-') tipoActual = IZQUIERDA;
        else if (alin2[k] == '-') tipoActual = ARRIBA;
        else tipoActual = DIAGONAL;

        if (tipoActual != DIAGONAL && tipoActual != tipoAnterior) {
            rupturas++;
        }
        tipoAnterior = tipoActual;
    }
    return rupturas;
}

// Dot plot
vector<vector<bool>> generarMatrizPuntos(const vector<char>& cadena1, const vector<char>& cadena2) {
    int n = (int)cadena1.size();
    int m = (int)cadena2.size();
    vector<vector<bool>> puntos(n, vector<bool>(m, false));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            puntos[i][j] = (cadena1[i] == cadena2[j]);
        }
    }
    return puntos;
}

void imprimirMatrizPuntosASCII(const vector<vector<bool>>& puntos, int maxFilas = 60, int maxColumnas = 100) {
    int n = (int)puntos.size();
    int m = n > 0 ? (int)puntos[0].size() : 0;
    if (n == 0 || m == 0) return;

    int filas = min(n, maxFilas);
    int columnas = min(m, maxColumnas);

    for (int i = 0; i < filas; i++) {
        string fila = "";
        for (int j = 0; j < columnas; j++) {
            fila += puntos[i][j] ? '*' : '.';
        }
        cout << fila << endl;
    }

    if (n > filas || m > columnas) {
        cout << "(" << filas << "x" << columnas << " de " << n << "x" << m << ")" << endl;
    }
}

map<string, vector<char>> leerSecuencias(const string& rutaArchivo) {
    map<string, vector<char>> secuencias;
    vector<string> nombresValidos = { "Bacteria", "Sars-Cov", "Influenza" };

    ifstream archivo(rutaArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo: " << rutaArchivo << endl;
        return secuencias;
    }

    string linea;
    string seccionActual = "";

    while (getline(archivo, linea)) {
        while (!linea.empty() && (linea.back() == '\r' || linea.back() == '\n')) {
            linea.pop_back();
        }

        size_t inicio = linea.find_first_not_of(" \t");
        if (inicio == string::npos) continue;
        size_t fin = linea.find_last_not_of(" \t");
        string trimmed = linea.substr(inicio, fin - inicio + 1);

        bool esEncabezado = false;
        for (auto& nombre : nombresValidos) {
            if (trimmed == nombre) {
                seccionActual = nombre;
                secuencias[seccionActual] = vector<char>();
                esEncabezado = true;
                break;
            }
        }
        if (esEncabezado) continue;
        if (seccionActual.empty()) continue;

        istringstream iss(trimmed);
        string token;
        while (iss >> token) {
            bool esNumero = !token.empty() &&
                all_of(token.begin(), token.end(), [](char c) { return isdigit((unsigned char)c); });
            if (esNumero) continue;

            for (char c : token) {
                if (isalpha((unsigned char)c)) {
                    secuencias[seccionActual].push_back((char)tolower((unsigned char)c));
                }
            }
        }
    }

    return secuencias;
}


void probarAlineamiento(const string& nombre1, vector<char>& cadena1,
    const string& nombre2, vector<char>& cadena2) {
    cout << "Alineando " << nombre1 << " (" << cadena1.size() << " nt) vs "
        << nombre2 << " (" << cadena2.size() << " nt)" << endl;

    auto inicio = chrono::high_resolution_clock::now();
    Matrix matriz = construccionMatriz(cadena1, cadena2);
    auto fin = chrono::high_resolution_clock::now();
    double segundos = chrono::duration<double>(fin - inicio).count();

    int scoreOptimo = matriz.Matriz[cadena1.size()][cadena2.size()].first;
    cout << "Score optimo: " << scoreOptimo << endl;
    cout << "Tiempo de construccion de la matriz: " << segundos << " s" << endl;

    size_t largoMuestra = min<size_t>(80, max(cadena1.size(), cadena2.size()) + 1);

    // Alineamiento
    pair<string, string> alineamientoIngenuo = reconstruirAlineamiento(matriz, cadena1, cadena2);
    int rupturasIngenuo = contarRupturas(alineamientoIngenuo.first, alineamientoIngenuo.second);
    cout << endl << "-- Primer camino encontrado --" << endl;
    cout << "Rupturas: " << rupturasIngenuo << endl;
    size_t muestraIngenuo = min(largoMuestra, alineamientoIngenuo.first.size());
    cout << alineamientoIngenuo.first.substr(0, muestraIngenuo) << endl;
    cout << alineamientoIngenuo.second.substr(0, muestraIngenuo) << endl;

    // Alineamiento con menos rupturas
    ResultadoAlineamiento mejor = seleccionarMenosRupturas(matriz, cadena1, cadena2);
    cout << endl << "-- Alineamiento con menos rupturas --" << endl;
    cout << "Rupturas: " << mejor.rupturas << endl;
    size_t muestraMejor = min(largoMuestra, mejor.alin1.size());
    cout << mejor.alin1.substr(0, muestraMejor) << endl;
    cout << mejor.alin2.substr(0, muestraMejor) << endl;

    // Numero total de alineamientos optimos
    vector<vector<long long>> memo(cadena1.size() + 1, vector<long long>(cadena2.size() + 1, -1));
    long long numCaminos = contarCaminos(matriz, cadena1.size(), cadena2.size(), memo);
    cout << endl;
    if (numCaminos == numeric_limits<long long>::max()) {
        cout << "Numero de alineamientos optimos: desborda long long" << endl;
    }
    else {
        cout << "Numero total de alineamientos optimos: " << numCaminos << endl;
    }

    // Dot plot
    vector<vector<bool>> puntos = generarMatrizPuntos(cadena1, cadena2);

    cout << endl << "-- Dot plot (filas=" << nombre1 << ", columnas=" << nombre2 << ") --" << endl;
    imprimirMatrizPuntosASCII(puntos);
    cout << endl;
}

int main() {
    map<string, vector<char>> secuencias = leerSecuencias("Sequencias.txt");

    if (secuencias.count("Bacteria") == 0 || secuencias.count("Sars-Cov") == 0 || secuencias.count("Influenza") == 0) {
        cerr << "No se pudieron leer correctamente las tres secuencias del archivo." << endl;
        return 1;
    }

    cout << "Longitudes leidas del archivo:" << endl;
    cout << " Bacteria : " << secuencias["Bacteria"].size() << " nucleotidos" << endl;
    cout << " Sars-Cov : " << secuencias["Sars-Cov"].size() << " nucleotidos" << endl;
    cout << " Influenza: " << secuencias["Influenza"].size() << " nucleotidos" << endl << endl;

    probarAlineamiento("Bacteria", secuencias["Bacteria"], "Sars-Cov", secuencias["Sars-Cov"]);
    probarAlineamiento("Bacteria", secuencias["Bacteria"], "Influenza", secuencias["Influenza"]);
    probarAlineamiento("Sars-Cov", secuencias["Sars-Cov"], "Influenza", secuencias["Influenza"]);

    return 0;
}