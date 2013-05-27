#include <iostream>
#include <Exception>
#include <cmath>
#include <stdlib.h>
#include "Amino.h"
#include <list>
#include <iomanip>
#include <windows.h>

void show(Amino*);
int fitness(Amino*);
void mutate(Amino*, double);
void sortHighscore(int, int*, Amino**);
long identifier(int, Amino*);


#pragma region Proteine
/*
const std::string PROT_1 = "1001001001001001001001";
const int FALTUNG_LENGTH = 22;
const int BEST_FITNESS = 99;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 0;
const double MUTATION_RATE_TOLERANCE = 4; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/

/*/*### PROTEIN 1 ######################################################
const std::string PROT_1 = "10100110100101100101";
const int FALTUNG_LENGTH = 20;
const int BEST_FITNESS = 9;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 1000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 50;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_TOLERANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/

/*### PROTEIN 2 ######################################################
const std::string PROT_1 = "0010011000011000011000011";
const int FALTUNG_LENGTH = 25;
const int BEST_FITNESS = 7;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 5;
const double MUTATION_RATE_TOLERANCE = 1; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/

//*### PROTEIN 3 ######################################################
const std::string PROT_1 = "000110011000001111111001100001100100";
const int FALTUNG_LENGTH = 36;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 13; // Meistens
//const int BEST_FITNESS = 14;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 200; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 8;
const double MUTATION_RATE_TOLERANCE = 8; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/

/*### PROTEIN 4 ######################################################
const std::string PROT_1 = "001001100110000011111111110000001100110010011111";
const int FALTUNG_LENGTH = 48;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 22;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 0;
const double MUTATION_RATE_TOLERANCE = 4; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/

/*### PROTEIN 5 ######################################################
const std::string PROT_1 = "11010101011110100010001000010001000101111010101011";
const int FALTUNG_LENGTH = 50;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 20;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 0;
const double MUTATION_RATE_TOLERANCE = 4; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf
//*/
#pragma endregion

// Anzahl der besten proteine, die unverändert in die nächste Generation übernommen werden
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;

// Elite Auswahl wird nur einmal übernommen, danach nur dessen Mutationen
const bool ELITE_ONLY_ONCE = true;

// Die Senkung der Mutationsrate pro Generation
const double MUTATION_RATE_LOWER = (MUTATION_RATE - MUTATION_RATE_LOWEST) / GENERATION_WITH_LOWEST_MUATTION_RATE;

// Anzeige
int NUMBER_WIDTH = 3;
const bool SHOW_DETAILS = false;
const bool SUCCESSES_ONLY = true;
const int MATRIX_SIZE = 20;

// Für Farbausgabe
HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );

void main()
{
	#pragma region Initialisierung und mehr
	time_t t;
    time(&t);
    srand((unsigned int)t);

	std::list<Amino*> populations = std::list<Amino*>();
	Amino* eliteProteins[ELITE_AMOUNT];

	// Ausgangsprotein
	Amino* proteinToCopy = new Amino[FALTUNG_LENGTH];

	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		proteinToCopy[i] = Amino(FORWARD, PROT_1[i] == 49);
	}

	// Titel
	SetConsoleTextAttribute(hstdout, 0x1f);
	std::cout << "Proteins -> ";
	if(SHOW_DETAILS)
	{
		for(int i = 0; i<AMOUNT_PROTEINS; i++) std::cout << std::setw(NUMBER_WIDTH) << std::right << i;
	}
	SetConsoleTextAttribute(hstdout, 0x0f);
	std::cout << " | Aver";
	std::cout << " Best";
	std::cout << " Best";
	std::cout << " Mut. Rate";
	std::cout << std::endl;
	#pragma endregion

	#pragma region Population generieren
	if(!SUCCESSES_ONLY) std::cout << std::setw(6) << std::right << "Gen    1 -> ";
	for(int i = 0; i<AMOUNT_PROTEINS; i++)
	{
		Amino* addMe = new Amino[FALTUNG_LENGTH];

		// Ausgangsprotein kopieren
		memcpy(addMe, proteinToCopy, FALTUNG_LENGTH * sizeof(*proteinToCopy));

		if(SHOW_DETAILS)
		{
			int f = fitness(addMe);
			
			std::cout << std::setw(NUMBER_WIDTH) << std::right << f;
		}

		populations.push_front(addMe);
	}

	std::cout << std::endl;
	#pragma endregion
	
	#pragma region In allen Generationen, alle Proteine durchgehen
	Amino* bestProtein = new Amino[FALTUNG_LENGTH];
	int bestProtein_fitness = -1;
	int bestProtein_round= -1;
	int scores[AMOUNT_PROTEINS];
	Amino* scores_values[AMOUNT_PROTEINS];
	for(int i = 0; i<GENERATIONS; i++)
	{
		bool showRow = !SUCCESSES_ONLY || i%500 == 0;

		if(showRow) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> ";

		#pragma region Mutieren
		double mutation_rate = MUTATION_RATE;

		// Toleranz in Mutationsrate einbauen
		double tolerance = 0;
		if(MUTATION_RATE_TOLERANCE != 0)
		{
			// Auf Mutation die Toleranz dazurechnen
			tolerance = (rand() % (int)(MUTATION_RATE_TOLERANCE * 10000) - (MUTATION_RATE_TOLERANCE/2.0 * 10000)) / 10000.0;

			mutation_rate += tolerance;
				
			// Auf Grenzen achten
			if(mutation_rate < MUTATION_RATE_LOWEST) mutation_rate = MUTATION_RATE_LOWEST;
		}

		// Mutieren lassen und Besten suchen
		std::list<Amino*>::iterator iter;
		double generationAverage = 0;
		int count = 0;
		for (iter = populations.begin(); iter != populations.end(); ++iter)
		{
			// Mutieren (aber nicht die "Elite", die bleibt bestehen)
			if(count >= ELITE_AMOUNT) mutate(*iter, mutation_rate);

			// Fitness ausrechnen
			int f = fitness(*iter);

			// Fitness in die Durschnittsrechnung einrechnen
			generationAverage += (f == -1 ? 0 : f);

			// Die Elite nicht in den Highscore (und damit in die nächste Generation übernehemen)
			if(count >= ELITE_AMOUNT || !ELITE_ONLY_ONCE)
			{
				// Fitness zu Highscore Array hinzufügen
				scores[count] = f;
				scores_values[count] = *iter;
			}

			// Fortlaufender index;
			count++;
		}

		// Highscore sortieren
		sortHighscore(AMOUNT_PROTEINS, scores, scores_values);

		// Auswerten, welches bestes Protein ist
		int bestGenerationProtein_fitness = scores[0];
		Amino* bestGenerationProtein = scores_values[0];

		// Durschnittsrechnung abschließen
		generationAverage = ceil(generationAverage / populations.size() * 10) / 10.0;
		#pragma endregion

		#pragma region Anzeigen
		if(SHOW_DETAILS)
		{
			for (iter = populations.begin(); iter != populations.end(); ++iter)
			{
				int f = fitness(*iter);

				if(f == -1) SetConsoleTextAttribute(hstdout, 0x0c);
				else if(*iter == bestGenerationProtein) SetConsoleTextAttribute(hstdout, 0xa0);

				if(showRow) std::cout << std::setw(NUMBER_WIDTH) << std::right << f;
			
				SetConsoleTextAttribute(hstdout, 0x0f);
			}
		}
		#pragma endregion
		
		#pragma region Auswerten, welches bestes Protein ist
		bool isNewBest = (bestGenerationProtein_fitness > bestProtein_fitness); 
		
		if(isNewBest)
		{
			memcpy(bestProtein, bestGenerationProtein, FALTUNG_LENGTH * sizeof(*bestGenerationProtein));
			bestProtein_fitness = bestGenerationProtein_fitness;
			bestProtein_round = i + 2;
		}

		// Auswertung anzeigen
		if(!showRow && isNewBest) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> "; // Wenn SUCCESSES_ONLY, dann hier einen Ersatz anzeigen

		if(showRow || isNewBest)
		{
			std::cout << " |" << std::setw(5) << generationAverage;
		
			std::cout << std::setw(5) << bestGenerationProtein_fitness;
		
			if(isNewBest) SetConsoleTextAttribute(hstdout, 0xa0);
			std::cout << std::setw(5) << bestProtein_fitness;
			SetConsoleTextAttribute(hstdout, 0x0f);

			std::cout << std::setw(10) << MUTATION_RATE;

			std::cout << std::endl;
		}
		#pragma endregion
		
		#pragma region Selektion
		// X Elite-Proteine wählen, die nicht gleich sind
		std::list<int> proteinIds = std::list<int>();
		for(int i = 0; i<AMOUNT_PROTEINS - ELITE_AMOUNT; i++)
		{
			// Id des Proteins bilden
			long id = identifier(FALTUNG_LENGTH, scores_values[i]);
			
			// Wenn Protein noch nicht in Liste vorhanden, dann hinzufügen
			std::list<int>::iterator findIter = std::find(proteinIds.begin(), proteinIds.end(), id);
			if(findIter == proteinIds.end())
			{
				eliteProteins[proteinIds.size()] = scores_values[i];
				proteinIds.push_back(id);
			}

			// Wenn Liste voll, dann abbrechen
			if(proteinIds.size() == ELITE_AMOUNT) break;
		}

		// Wenn zu wenig Elite Proteine vorhanden, dann Rest mit Bestem Protein auffüllen
		if(proteinIds.size() < ELITE_AMOUNT)
		{
			for(int i = proteinIds.size(); i<ELITE_AMOUNT; i++) eliteProteins[i] = scores_values[0];
		}

		// Neue Population generieren
		std::list<Amino*> newPopulation = std::list<Amino*>();
		for(int i = 0; i<AMOUNT_PROTEINS; i++)
		{
			// Wenn es eine Elite gibt, dann aus dieser die zu Übernehmenden wählen
			int copyFrom = (ELITE_AMOUNT != 0 ? (i % ELITE_AMOUNT) : i);

			Amino* newProt = new Amino[FALTUNG_LENGTH];
			memcpy(newProt, eliteProteins[copyFrom], FALTUNG_LENGTH * sizeof(*newProt));

			newPopulation.push_back(newProt);
		}

		populations = newPopulation;
		#pragma endregion

		// Wenn Maximum erreicht, dann abbrechen
		if(bestProtein_fitness == BEST_FITNESS) break;

		// Mutationsrate senken
		MUTATION_RATE -= MUTATION_RATE_LOWER;

		if(MUTATION_RATE < MUTATION_RATE_LOWEST) MUTATION_RATE = MUTATION_RATE_LOWEST;
	}
	
	std::cout << std::endl;
	#pragma endregion

	#pragma region Ausgabe des Highscores
	SetConsoleTextAttribute(hstdout, 0x70);
	std::cout << std::setw(80) << std::left << " Highscore" << std::endl;
	SetConsoleTextAttribute(hstdout, 0x0f);
	std::cout << "Best protein reached in generation " << bestProtein_round << std::endl;
	for(int i = 0; i<ELITE_AMOUNT && i<5; i++)
	{
		int f = fitness(eliteProteins[i]);
		int id = identifier(FALTUNG_LENGTH, eliteProteins[i]);
		std::cout << "Rank " << (i+1) << " (Fitness: " << f << ", Id: " << id << ")" << std::endl;
		show(eliteProteins[i]);
	}
	#pragma endregion
}

void mutate(Amino* faltung, double mutationRate)
{
	// Anzahl der Mutationen berechnen
	int mutations = (mutationRate / 100.0) * FALTUNG_LENGTH;

	// Verhindern, dass Mutationsrate so klein ist, dass gar keine Mutataion mehr stattfindet
	mutations = (mutations < 1 ? 1 : mutations);

	for(int i=0; i<mutations; i++)
	{
		// Nach Zufall eine Aminosäure ermitteln
		int randAmino = rand() % FALTUNG_LENGTH;
	
		// so lange eine neue zufällige Richtung ermitteln, bis eine NEUE Richtung ermittelt wurde
		int randDir;
		do
		{
			randDir = rand() % 3;

			// Zufallszahl in die Richtung aus dem enum umwandeln
			switch(randDir)
			{
				case 0: randDir = FORWARD; break;
				case 1: randDir = LEFT; break;
				case 2: randDir = RIGHT; break;
			}
		} while(randDir == faltung[randAmino].GetDir());

		// Richtung in der Aminosäure setzen
		faltung[randAmino].SetDir((Direction)randDir);
	}
}

// Generates a number, that identifies a protein
long identifier(int length, Amino* protein)
{
	int id = 0;

	for(int i = 0; i<length; i++)
	{
		int aminoId = 0;
		
		if(protein[i].GetDir() == LEFT) aminoId = 0;
		else if(protein[i].GetDir() == FORWARD) aminoId = 1;
		else if(protein[i].GetDir() == RIGHT) aminoId = 2;
		
		if(protein[i].IsHydrophob()) aminoId += 3;

		//*
		id += aminoId * (i+1) * 100;
		/*/
		id = id << 3;
		id += aminoId;
		//*/
	}

	return id;
}

void sortHighscore(int length, int* scores, Amino** score_values)
{
	int temp;
	Amino* temp_amino;
   
	for(int i=0; i<length; i++)
	{
		for(int j=0; j<length-1; j++)
		{
			if(scores[j] < scores[j+1])
			{
				temp = scores[j];
				scores[j] = scores[j+1];
				scores[j+1] = temp;

				temp_amino = score_values[j];
				score_values[j] = score_values[j+1];
				score_values[j+1] = temp_amino;
			}
		}
	}
}

int fitness(Amino* faltung)
{
	int x, y;
	const int width = MATRIX_SIZE;
	const int height = MATRIX_SIZE;
	int matrix[width][height];
	int startX = width / 2, startY = height / 2;

	int connections = 0;
	int directConnections = 0;
	bool lastIsHydrophob = false;

	try
	{
		#pragma region Initialisierung
		// Alle Felder mit -1 (==leer) initialisieren
		for(int x = 0; x<height; x++)
		{
			for(int y = 0; y<width; y++)
			{
				matrix[x][y] = -1;
			}
		}

		// Anfangskoordinaten setzen
		x = startX;
		y = startY;
		#pragma endregion

		#pragma region Aminos in Matrix eintragen
		int dir = NORTH;
		int maxX = startX, minX = startX, maxY = startY, minY= startY;
		for(int i = 0; i<FALTUNG_LENGTH; i++)
		{
			Amino* amino = &faltung[i];

			// Die richtige Himmelsrichtung einschlagen
			switch(amino->GetDir())
			{
				// Bei Vorwärts ändert sich nichts
				case FORWARD: break;
				// Ansonsten ändert sich die Richtung
				case LEFT: dir = (dir == 0 ? 3: (dir-1)); break;
				case RIGHT: dir = (dir + 1) % 4; break;
			}

			// Prüfen, ob sich Aminos überlappen
			switch(dir)
			{
				case NORTH:
					y++;
					break;
				case SOUTH:
					y--;
					break;
				case EAST:
					x++;
					break;
				case WEST:
					x--;
					break;
			}
			//std::cout << x << "/" << y << std::endl;
			if(matrix[x][y] != -1) throw -1;

			// Grenzen ermitteln
			if(x > maxX) maxX = x;
			if(x < minX) minX = x;
			if(y > maxY) maxY = y;
			if(y < minY) minY = y;

			// Direkte Verbindungen zählen
			if(lastIsHydrophob && amino->IsHydrophob()) directConnections++;

			// In Matrix eintragen
			matrix[x][y] = (amino->IsHydrophob() ? 1 : 0);

			lastIsHydrophob = amino->IsHydrophob();
		}
		#pragma endregion

		#pragma region Alle Verbindungen in der Matrix zählen
		for(int x = 0; x<height; x++)
		{
			for(int y = 0; y<width; y++)
			{
				int hydrophob = matrix[x][y];

				// Wenn leer, dann weiter
				if(hydrophob == -1) continue;
				// Wenn Hydrophob, dann beahndeln
				else if(hydrophob == 1)
				{
					// Alle Verbindungen des Hydrophobs zählen
					if(matrix[x+1][y] == 1) connections++;
					if(matrix[x-1][y] == 1) connections++;
					if(matrix[x][y+1] == 1) connections++;
					if(matrix[x][y-1] == 1) connections++;

					// Markierung auf leer setzen
					matrix[x][y] = -1;
				}
			}
		}
		#pragma endregion
	}
	catch(int e)
	{
		return e;
	}

	// Alle drekten Verbundungen von der Gesamtzahl aller Verbindungen abziehen
	int f = connections - directConnections;

	// Wenn die gezählte Zahl kleiner 0 ist, dann liegt ein Fehler vor, der mit einer -2 signalisiert wird
	if(f < 0) return -2;

	// Fitness zurückgeben
	return f;
}

void show(Amino* faltung)
{
	int x, y;
	const int width = 200;
	const int height = 200;
	Amino::DisplayInfo* matrix[width][height];
	int startX = width / 2, startY = height / 2;

	#pragma region Initialisierung
	// Alle Felder mit -1 (==leer) initialisieren
	for(int x = 0; x<height; x++)
	{
		for(int y = 0; y<width; y++)
		{
			matrix[x][y] = nullptr;
		}
	}
	//*/

	// Anfangskoordinaten setzen
	x = startX;
	y = startY;
	#pragma endregion

	#pragma region Aminos in Matrix eintragen
	Amino::DisplayInfo* lastAminoDisplayInfo = nullptr;
	int dirNr = NORTH;
	SkyDirection sDir = NORTH;
	int maxX = startX, minX = startX, maxY = startY, minY= startY;
	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		Amino* amino = &faltung[i];

		// Die richtige Himmelsrichtung einschlagen
		switch(amino->GetDir())
		{
			// Bei Vorwärts ändert sich nichts
			case FORWARD: break;
			// Ansonsten ändert sich die Richtung
			case LEFT: dirNr = (sDir == 0 ? 3: (sDir-1)); break;
			case RIGHT: dirNr = (sDir + 1) % 4; break;
		}

		// std::cout << (amino->GetDir() == FORWARD ? "FORWARD" : (amino->GetDir() == LEFT ? "LEFT" : "RIGHT")) << " | ";
			
		switch(dirNr)
		{
			case NORTH:
				sDir = NORTH;
				y++;
				break;
			case SOUTH:
				sDir = SOUTH;
				y--;
				break;
			case EAST:
				sDir = EAST;
				x++;
				break;
			case WEST:
				sDir = WEST;
				x--;
				break;
		}
		
		// Grenzen ermitteln
		if(x > maxX) maxX = x;
		if(x < minX) minX = x;
		if(y > maxY) maxY = y;
		if(y < minY) minY = y;

		if(lastAminoDisplayInfo != nullptr) lastAminoDisplayInfo->SetRelDir(amino->GetDir());

		lastAminoDisplayInfo = new Amino::DisplayInfo(sDir, amino->IsHydrophob());
		
		// In Matrix eintragen
		matrix[x][y] = lastAminoDisplayInfo;
	}
	#pragma endregion
	
	#pragma region Aminos anzeigen	
	//const int protWidth = maxX - minX;
	//const int protHeight = maxY - minY;

	maxX++;
	maxY++;
	minX--;
	minY--;

	//std::cout << "Protein: (X: " << minX << "/" << maxX << " | Y: " << minY << "/" << maxY << ")" << std::endl;
	
	for(int y = maxY; y>=minY; y--)
	{
		std::cout << std::setw(3) << (-(y-maxY)+1);

		for(int x = minX; x<=maxX; x++)
		{
			Amino::DisplayInfo* amino = matrix[x][y];

			// Wenn leer, dann weiter
			if(amino == nullptr) std::cout << ' ';
			else
			{
				if(amino->IsHydrophob()) SetConsoleTextAttribute(hstdout, (x%2 + y%2 == 1 ? 0x2f : 0xaf));
				else SetConsoleTextAttribute(hstdout, 0x0f);

				Direction relDir = amino->GetRelDir();

				switch(amino->GetSkyDir())
				{
					case NORTH:
						if(relDir == FORWARD) std::cout << (char)179;
						else if(relDir == LEFT) std::cout << (char)191;
						else std::cout << (char)218;
						break;
					case SOUTH:
						if(relDir == FORWARD) std::cout << (char)179;
						else if(relDir == LEFT) std::cout << (char)192;
						else std::cout << (char)217;
						break;
					case EAST:
						if(relDir == FORWARD) std::cout << (char)196;
						else if(relDir == LEFT) std::cout << (char)217;
						else std::cout << (char)191;
						break;
					case WEST:
						if(relDir == FORWARD) std::cout << (char)196;
						else if(relDir == LEFT) std::cout << (char)218;
						else std::cout << (char)192;
						break;
				}

				SetConsoleTextAttribute(hstdout, 0x0f);
			}
		}

		std::cout << std::endl;
	}
	#pragma endregion
}