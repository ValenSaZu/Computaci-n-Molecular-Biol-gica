#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <optional>
#include <functional>
#include <cctype>

using namespace std;

struct EsquemaScore {
    function<int(char, char)> sustitucion = [](char a, char b) {
        return a == b ? 1 : -1;
    };
    int gap = -2;
};

struct ResultadoPar {
    string alin1, alin2;
    int score = 0;
};

ResultadoPar alinearPar(string& s1, string& s2, EsquemaScore& esquema = EsquemaScore{}) {
    size_t n = s1.size(), m = s2.size();
    vector<vector<int>> M(n + 1, vector<int>(m + 1, 0));

    for (size_t i = 1; i <= n; ++i) M[i][0] = static_cast<int>(i) * esquema.gap;
    for (size_t j = 1; j <= m; ++j) M[0][j] = static_cast<int>(j) * esquema.gap;

    for (size_t i = 1; i <= n; ++i)
        for (size_t j = 1; j <= m; ++j) {
            int diag   = M[i - 1][j - 1] + esquema.sustitucion(s1[i - 1], s2[j - 1]);
            int arriba = M[i - 1][j] + esquema.gap;
            int izq    = M[i][j - 1] + esquema.gap;
            M[i][j] = max({diag, arriba, izq});
        }

    string a1, a2;
    size_t i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && M[i][j] == M[i - 1][j - 1] + esquema.sustitucion(s1[i - 1], s2[j - 1])) {
            a1 = s1[i - 1] + a1; a2 = s2[j - 1] + a2; --i; --j;
        } else if (i > 0 && M[i][j] == M[i - 1][j] + esquema.gap) {
            a1 = s1[i - 1] + a1; a2 = '-' + a2; --i;
        } else {
            a1 = '-' + a1; a2 = s2[j - 1] + a2; --j;
        }
    }
    return {a1, a2, M[n][m]};
}

map<string, string> leerPrimersBRCA1(string& ruta) {
    map<string, string> secuencias;

    ifstream archivo(ruta);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo: " << ruta << '\n';
        return secuencias;
    }

    string linea, regionActual;
    while (getline(archivo, linea)) {
        while (!linea.empty() && (linea.back() == '\r' || linea.back() == '\n' || linea.back() == '\t'))
            linea.pop_back();

        size_t ini = linea.find_first_not_of(" \t");
        if (ini == string::npos) continue;
        string trimmed = linea.substr(ini);

        istringstream iss(trimmed);
        string primerToken, segundoToken;
        iss >> primerToken;
        if (primerToken.size() == 1 && isupper((unsigned char)primerToken[0])) {
            iss >> segundoToken;
            if (segundoToken == "BRCA1") regionActual = primerToken;
        }

        char tipo = 0;
        if (trimmed.find("F:") != string::npos) tipo = 'F';
        else if (trimmed.find("R:") != string::npos) tipo = 'R';
        if (tipo == 0 || regionActual.empty()) continue;

        string mejor, actual;
        for (char c : trimmed) {
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') {
                actual += c;
            } else {
                if (actual.size() > mejor.size()) mejor = actual;
                actual.clear();
            }
        }
        if (actual.size() > mejor.size()) mejor = actual;
        if (mejor.size() < 8) continue;

        secuencias[regionActual + "_" + tipo] = mejor;
    }
    return secuencias;
}

vector<vector<int>> matrizDeScores(vector<string>& secs, EsquemaScore& esquema) {
    size_t k = secs.size();
    vector<vector<int>> S(k, vector<int>(k, 0));
    for (size_t i = 0; i < k; ++i)
        for (size_t j = i + 1; j < k; ++j) {
            auto [a1, a2, score] = alinearPar(secs[i], secs[j], esquema);
            S[i][j] = S[j][i] = score;
        }
    return S;
}

size_t encontrarCentro(vector<vector<int>>& S) {
    vector<int> sumas(S.size());
    for (size_t i = 0; i < S.size(); ++i) sumas[i] = accumulate(S[i].begin(), S[i].end(), 0);

    auto it = max_element(sumas.begin(), sumas.end());
    return static_cast<size_t>(distance(sumas.begin(), it));
}

vector<int> insercionesRelativas(string& centroOriginal, string& centroAlineado) {
    vector<int> ins(centroOriginal.size() + 1, 0);
    size_t j = 0;
    for (char c : centroAlineado) {
        if (c == '-') ++ins[j];
        else ++j;
    }
    return ins;
}

vector<int> combinarInserciones(vector<vector<int>>& todas) {
    vector<int> global = todas.front();
    for (auto& v : todas) transform(global.begin(), global.end(), v.begin(), global.begin(), [](int a, int b) { return max(a, b); });
    return global;
}

string reconstruirCentro(string& centro, vector<int>& globalIns) {
    string fila;
    for (size_t j = 0; j < centro.size(); ++j) {
        fila.append(globalIns[j], '-');
        fila += centro[j];
    }
    fila.append(globalIns[centro.size()], '-');
    return fila;
}

string reconstruirFila(string& centroAlineado, string& otroAlineado,
                        vector<int>& globalIns) {
    string fila;
    size_t j = 0, usados = 0;
    for (size_t p = 0; p < centroAlineado.size(); ++p) {
        if (centroAlineado[p] == '-') {
            fila += otroAlineado[p];
            ++usados;
        } else {
            fila.append(globalIns[j] - usados, '-');
            fila += otroAlineado[p];
            ++j;
            usados = 0;
        }
    }
    fila.append(globalIns[j] - usados, '-');
    return fila;
}

struct AlineamientoMultiple {
    vector<string> nombres;
    vector<string> filas;
};

AlineamientoMultiple alineamientoEstrella(map<string, string>& secuenciasPorNombre, EsquemaScore& esquema = EsquemaScore{}) {
    vector<string> nombres, secs;
    for (auto& [nombre, seq] : secuenciasPorNombre) {
        nombres.push_back(nombre);
        secs.push_back(seq);
    }

    size_t k = secs.size();
    AlineamientoMultiple resultado{nombres, secs};
    if (k < 2) return resultado;

    auto S = matrizDeScores(secs, esquema);
    size_t c = encontrarCentro(S);

    vector<optional<ResultadoPar>> pares(k);
    for (size_t i = 0; i < k; ++i) if (i != c) pares[i] = alinearPar(secs[c], secs[i], esquema);

    vector<vector<int>> todasInserciones;
    for (size_t i = 0; i < k; ++i) if (i != c) todasInserciones.push_back(insercionesRelativas(secs[c], pares[i]->alin1));
    vector<int> globalIns = combinarInserciones(todasInserciones);

    vector<string> filas(k);
    filas[c] = reconstruirCentro(secs[c], globalIns);
    for (size_t i = 0; i < k; ++i) if (i != c) filas[i] = reconstruirFila(pares[i]->alin1, pares[i]->alin2, globalIns);

    cout << "Secuencia centro elegida: " << nombres[c] << " (suma de scores = " << accumulate(S[c].begin(), S[c].end(), 0) << ")\n\n";

    return {nombres, filas};
}

void imprimirAlineamiento(AlineamientoMultiple& msa) {
    size_t ancho = 0;
    for (auto& n : msa.nombres) ancho = max(ancho, n.size());
    for (size_t i = 0; i < msa.nombres.size(); ++i) cout << msa.nombres[i] << string(ancho - msa.nombres[i].size() + 1, ' ') << msa.filas[i] << '\n';
}

int main() {
    map<string, string> ejemploDiapositiva = {
        {"S1", "ATTGCCATT"},
        {"S2", "ATGGCCATT"},
        {"S3", "ATCCAATTTT"},
        {"S4", "ATCTTCTT"},
        {"S5", "ACTGACC"}
    };

    cout << "Ejemplo de las diapositivas (S1..S5)\n";
    cout << "Numero de secuencias a alinear (k): " << ejemploDiapositiva.size() << "\n\n";
    imprimirAlineamiento(alineamientoEstrella(ejemploDiapositiva));

    cout << "\nPrimers BRCA1 (BRCA1.txt) - todos juntos\n";
    map<string, string> primersBRCA1 = leerPrimersBRCA1("BRCA1.txt");
    if (primersBRCA1.size() < 2) {
        cerr << "No se pudieron leer al menos 2 secuencias de BRCA1.txt\n";
        return 1;
    }
    cout << "Numero de secuencias a alinear (k): " << primersBRCA1.size() << "\n\n";
    imprimirAlineamiento(alineamientoEstrella(primersBRCA1));

    map<string, string> soloForward, soloReverse;
    for (auto& [nombre, seq] : primersBRCA1) {
        if (nombre.size() >= 2 && nombre.back() == 'F') soloForward[nombre] = seq;
        else if (nombre.size() >= 2 && nombre.back() == 'R') soloReverse[nombre] = seq;
    }

    cout << "\n Solo cadenas Forward (" << soloForward.size() << ") \n";
    imprimirAlineamiento(alineamientoEstrella(soloForward));

    cout << "\n Solo cadenas Reverse (" << soloReverse.size() << ") \n";
    imprimirAlineamiento(alineamientoEstrella(soloReverse));

    return 0;
}