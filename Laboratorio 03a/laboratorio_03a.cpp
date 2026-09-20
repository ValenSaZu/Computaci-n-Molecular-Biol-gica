#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <climits>
#include <cmath>

using namespace std;

char complemento(char base) {
    switch (base) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default:  return base;
    }
}

string complementoReverso(const string& cadena) {
    string resultado(cadena.rbegin(), cadena.rend());
    for (char& c : resultado) c = complemento(c);
    return resultado;
}

int overlap(const string& a, const string& b) {
    int maxT = (int)min(a.size(), b.size());
    for (int t = maxT; t >= 0; t--) {
        if (a.compare(a.size() - t, t, b, 0, t) == 0) {
            return t;
        }
    }
    return 0;
}

struct Fragmento {
    string nombre;
    string directa;
    string reversa;
};

struct NodoOrientado {
    int idFragmento;
    int orientacion;
    string secuencia;
};

vector<NodoOrientado> construirNodosOrientados(const vector<Fragmento>& frags) {
    vector<NodoOrientado> nodos;
    for (int i = 0; i < (int)frags.size(); i++) {
        nodos.push_back({ i, 0, frags[i].directa });
        nodos.push_back({ i, 1, frags[i].reversa });
    }
    return nodos;
}

struct Ensamble {
    vector<int> ordenIds;
    vector<int> orientaciones;
    vector<int> pesosAristas;
    int pesoTotal = -1;
    int longitudConsenso = -1;
};

struct BuscadorConsenso {
    const vector<Fragmento>& frags;
    vector<NodoOrientado> nodos;
    vector<vector<int>> pesoArista;
    int n;
    long long sumaLongitudes;
    int longitudObjetivo;

    Ensamble mejor;
    int mejorDiferencia = INT_MAX;

    vector<int> ordenActual, orientActual, pesosActual;
    long long pesoActual = 0;

    BuscadorConsenso(const vector<Fragmento>& f, int objetivo)
        : frags(f), longitudObjetivo(objetivo) {
        n = (int)frags.size();
        nodos = construirNodosOrientados(frags);
        int total = (int)nodos.size();
        pesoArista.assign(total, vector<int>(total, -1));
        for (int u = 0; u < total; u++) {
            for (int v = 0; v < total; v++) {
                if (nodos[u].idFragmento == nodos[v].idFragmento) continue;
                pesoArista[u][v] = overlap(nodos[u].secuencia, nodos[v].secuencia);
            }
        }
        sumaLongitudes = 0;
        for (auto& f2 : frags) sumaLongitudes += (long long)f2.directa.size();
    }

    void evaluarHoja() {
        int longitudConsenso = (int)(sumaLongitudes - pesoActual);
        int diferencia = abs(longitudConsenso - longitudObjetivo);
        if (diferencia < mejorDiferencia) {
            mejorDiferencia = diferencia;
            mejor.ordenIds = ordenActual;
            mejor.orientaciones = orientActual;
            mejor.pesosAristas = pesosActual;
            mejor.pesoTotal = (int)pesoActual;
            mejor.longitudConsenso = longitudConsenso;
        }
    }

    void backtrack(vector<bool>& usado, int ultimoNodo) {
        if ((int)ordenActual.size() == n) {
            evaluarHoja();
            return;
        }
        for (int id = 0; id < n; id++) {
            if (usado[id]) continue;
            for (int orient = 0; orient < 2; orient++) {
                int nodoActual = id * 2 + orient;
                int pesoArco = 0;
                if (!ordenActual.empty()) {
                    pesoArco = pesoArista[ultimoNodo][nodoActual];
                }
                bool esPrimerNodo = ordenActual.empty();
                usado[id] = true;
                ordenActual.push_back(id);
                orientActual.push_back(orient);
                if (!esPrimerNodo) pesosActual.push_back(pesoArco);
                pesoActual += pesoArco;

                backtrack(usado, nodoActual);

                pesoActual -= pesoArco;
                if (!esPrimerNodo) pesosActual.pop_back();
                ordenActual.pop_back();
                orientActual.pop_back();
                usado[id] = false;
            }
        }
    }

    Ensamble buscar() {
        vector<bool> usado(n, false);
        backtrack(usado, -1);
        return mejor;
    }
};

struct Layout {
    vector<string> filas;
    vector<string> nombres;
    vector<int> inicio;
    string consenso;
};

Layout construirLayout(const vector<Fragmento>& frags, const Ensamble& ens) {
    int n = (int)ens.ordenIds.size();
    vector<string> secuenciasUsadas(n);
    vector<int> inicio(n);

    inicio[0] = 0;
    for (int k = 0; k < n; k++) {
        int id = ens.ordenIds[k];
        secuenciasUsadas[k] = (ens.orientaciones[k] == 0) ? frags[id].directa
            : frags[id].reversa;
    }
    for (int k = 1; k < n; k++) {
        int solape = ens.pesosAristas[k - 1];
        inicio[k] = inicio[k - 1] + (int)secuenciasUsadas[k - 1].size() - solape;
    }

    int longitudTotal = inicio[n - 1] + (int)secuenciasUsadas[n - 1].size();

    string consenso(longitudTotal, '?');
    for (int k = 0; k < n; k++) {
        for (int p = 0; p < (int)secuenciasUsadas[k].size(); p++) {
            int col = inicio[k] + p;
            char base = secuenciasUsadas[k][p];
            if (consenso[col] == '?') consenso[col] = base;
            else if (consenso[col] != base) {
                cerr << "Aviso: discrepancia en columna " << col
                    << " (" << consenso[col] << " vs " << base << ")" << endl;
            }
        }
    }

    Layout layout;
    layout.consenso = consenso;
    layout.inicio = inicio;
    for (int k = 0; k < n; k++) {
        string fila(longitudTotal, '-');
        for (int p = 0; p < (int)secuenciasUsadas[k].size(); p++) {
            fila[inicio[k] + p] = secuenciasUsadas[k][p];
        }
        layout.filas.push_back(fila);
        int id = ens.ordenIds[k];
        string etiqueta = frags[id].nombre + (ens.orientaciones[k] == 1 ? "'" : "");
        layout.nombres.push_back(etiqueta);
    }
    return layout;
}

void imprimirLayout(const Layout& layout) {
    for (size_t k = 0; k < layout.filas.size(); k++) {
        cout << layout.nombres[k] << ": " << layout.filas[k] << endl;
    }
    cout << string(layout.nombres[0].size() + 2, ' ') << layout.consenso << "  <- consenso" << endl;
}

struct GrafoFiltrado {
    int n;
    vector<vector<int>> pesoUsable;
};

GrafoFiltrado construirGrafoFiltrado(const vector<string>& secuenciasOrientadas, int t) {
    int n = (int)secuenciasOrientadas.size();
    GrafoFiltrado g;
    g.n = n;
    g.pesoUsable.assign(n, vector<int>(n, -1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            int w = overlap(secuenciasOrientadas[i], secuenciasOrientadas[j]);
            if (w >= t) g.pesoUsable[i][j] = w;
        }
    }
    return g;
}

void buscarHamiltonianos(const GrafoFiltrado& g, vector<int>& camino,
    vector<bool>& usado, vector<vector<int>>& resultados) {
    if ((int)camino.size() == g.n) {
        resultados.push_back(camino);
        return;
    }
    for (int v = 0; v < g.n; v++) {
        if (usado[v]) continue;
        if (!camino.empty() && g.pesoUsable[camino.back()][v] < 0) continue;
        usado[v] = true;
        camino.push_back(v);
        buscarHamiltonianos(g, camino, usado, resultados);
        camino.pop_back();
        usado[v] = false;
    }
}

struct UnionFind {
    vector<int> padre;
    UnionFind(int n) : padre(n) { for (int i = 0; i < n; i++) padre[i] = i; }
    int find(int x) { return padre[x] == x ? x : padre[x] = find(padre[x]); }
    void unir(int a, int b) { padre[find(a)] = find(b); }
};

void reportarParticionMulticontig(const GrafoFiltrado& g, const vector<string>& etiquetas) {
    UnionFind uf(g.n);
    for (int i = 0; i < g.n; i++)
        for (int j = 0; j < g.n; j++)
            if (g.pesoUsable[i][j] >= 0) uf.unir(i, j);

    vector<vector<int>> componentes(g.n);
    for (int i = 0; i < g.n; i++) componentes[uf.find(i)].push_back(i);

    int numeroContig = 0;
    for (int raiz = 0; raiz < g.n; raiz++) {
        if (componentes[raiz].empty()) continue;
        numeroContig++;
        vector<int>& miembros = componentes[raiz];

        cout << "  Contig " << numeroContig << " (" << miembros.size() << " fragmentos): ";
        for (int idx : miembros) cout << etiquetas[idx] << " ";
        cout << endl;

        if (miembros.size() == 1) continue;

        int m = (int)miembros.size();
        GrafoFiltrado sub;
        sub.n = m;
        sub.pesoUsable.assign(m, vector<int>(m, -1));
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                sub.pesoUsable[a][b] = g.pesoUsable[miembros[a]][miembros[b]];

        vector<int> camino;
        vector<bool> usado(m, false);
        vector<vector<int>> resultados;
        buscarHamiltonianos(sub, camino, usado, resultados);

        cout << "    Caminos Hamiltonianos dentro del contig: " << resultados.size() << endl;
        for (auto& r : resultados) {
            cout << "      ";
            for (int idx : r) cout << etiquetas[miembros[idx]] << " ";
            cout << endl;
        }
    }
}

int main() {
    vector<Fragmento> frags = {
        {"f1", "ATCCGTTGAAGCCGCGGGC", ""},
        {"f2", "TTAACTCGAGG", ""},
        {"f3", "TTAAGTACTGCCCG", ""},
        {"f4", "ATCTGTGTCGGG", ""},
        {"f5", "CGACTCCCGACACA", ""},
        {"f6", "CACAGATCCGTTGAAGCCGCGGG", ""},
        {"f7", "CTCGAGTTAAGTA", ""},
        {"f8", "CGCGGGCAGTACTT", ""}
    };
    for (auto& f : frags) f.reversa = complementoReverso(f.directa);

    cout << " 2.1 SECUENCIA DE CONSENSO" << endl;

    long long sumaLongitudes = 0;
    for (auto& f : frags) {
        cout << f.nombre << " (directa) : " << f.directa
            << "  [len=" << f.directa.size() << "]" << endl;
        cout << f.nombre << " (reversa) : " << f.reversa << endl;
        sumaLongitudes += (long long)f.directa.size();
    }
    cout << "Suma de longitudes ||F|| = " << sumaLongitudes << endl;

    int longitudObjetivo = 55;
    cout << "Longitud objetivo l = " << longitudObjetivo << endl << endl;

    BuscadorConsenso buscador(frags, longitudObjetivo);
    Ensamble mejor = buscador.buscar();

    cout << "Mejor ensamble encontrado:" << endl;
    cout << "  Orden       : ";
    for (int k = 0; k < (int)mejor.ordenIds.size(); k++)
        cout << frags[mejor.ordenIds[k]].nombre
        << (mejor.orientaciones[k] == 1 ? "'" : "") << " ";
    cout << endl;
    cout << "  Pesos aristas: ";
    for (int w : mejor.pesosAristas) cout << w << " ";
    cout << endl;
    cout << "  Peso total w(P)     = " << mejor.pesoTotal << endl;
    cout << "  Longitud consenso   = " << mejor.longitudConsenso << endl;
    cout << "  |diferencia con l|  = " << abs(mejor.longitudConsenso - longitudObjetivo) << endl;
    cout << endl;

    Layout layout = construirLayout(frags, mejor);
    cout << "Alineamiento / layout:" << endl;
    imprimirLayout(layout);
    cout << endl;
    cout << "Secuencia de consenso final (longitud " << layout.consenso.size() << "):" << endl;
    cout << layout.consenso << endl;

    int linkageMinimo = *min_element(mejor.pesosAristas.begin(), mejor.pesosAristas.end());
    cout << endl << "Link mas debil del ensamble (linkage minimo observado) = "
        << linkageMinimo << endl;

    cout << " 2.2 SUBGRAFOS ACICLICOS (busqueda con linkage t)" << endl;

    int n = (int)mejor.ordenIds.size();
    vector<string> secOrientadas(n);
    vector<string> etiquetas(n);
    for (int k = 0; k < n; k++) {
        int id = mejor.ordenIds[k];
        secOrientadas[k] = (mejor.orientaciones[k] == 0) ? frags[id].directa : frags[id].reversa;
        etiquetas[k] = frags[id].nombre + (mejor.orientaciones[k] == 1 ? "'" : "");
    }

    vector<int> valoresT = { linkageMinimo, linkageMinimo + 1, 5, 8, 11, 13, 19 };
    for (int t : valoresT) {
        cout << endl << "--- t = " << t << " ---" << endl;
        GrafoFiltrado g = construirGrafoFiltrado(secOrientadas, t);

        cout << "Aristas del grafo OM(F," << t << "):" << endl;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (g.pesoUsable[i][j] >= 0)
                    cout << "  " << etiquetas[i] << " -> " << etiquetas[j]
                    << "  (peso " << g.pesoUsable[i][j] << ")" << endl;

        vector<int> camino;
        vector<bool> usado(n, false);
        vector<vector<int>> resultados;
        buscarHamiltonianos(g, camino, usado, resultados);

        cout << "Caminos Hamiltonianos encontrados: " << resultados.size() << endl;
        for (size_t r = 0; r < resultados.size(); r++) {
            cout << "  Camino " << (r + 1) << ": ";
            for (int idx : resultados[r]) cout << etiquetas[idx] << " ";
            cout << endl;
        }
        if (resultados.empty()) {
            cout << "No hay un unico contig: particion (Multicontig) resultante:" << endl;
            reportarParticionMulticontig(g, etiquetas);
        }
    }

    return 0;
}