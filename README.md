# TANSProject
Progetto per l'esame di Tecniche di Analisi Numerica e Simulazione.

Per usare il programma:
1) Aprire una sessione di ROOT;
2) Compilare il file "LibraryCompiler.C"
3) Eseguire la funzione "LibraryCompiler()"
4) Eseguire la funzione "Simulation()"
E' definita come

Simulation(int N_exp = 1e6, unsigned int seed = 69420, int zgen_flag = 1, int multiplicity_flag = 1, int multiscattering_flag = 1, int smearing_flag = 1, int N_noise = 0, const char* input_file = "kinem.root", const char* output_file = "simulation.root")

Valori di multiplicity_flag: 1 per la distribuzione data, 2 per un valore costante, 3 per distribuzione uniforme
Valori di zgen_flag: 1 per distribuzione gaussiana, 2 per distribuzione uniforme

5) Eseguire la funzione "Reconstruction()"

Per fare simulazioni o ricostruzioni non standard, seguire le istruzioni inserite come commenti nei file "LampadinaSimulation.C" e "LampadinaReconstrucion.C"
