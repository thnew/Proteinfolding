#include <iostream>
#include <Exception>
#include <cmath>
#include <stdlib.h>
#include "Protein.h"
#include "Amino.h"
#include <list>
#include <iomanip>
#include <string>
#include <windows.h>

// für Quicksort
int compare (const void * a, const void * b)
{
	Protein* protA = *(Protein**)a;
	Protein* protB = *(Protein**)b;

	int returnMe = (protB->CalcFitness() - protA->CalcFitness());

	return returnMe;
}

#pragma region Proteine
/*### Testprotein ######################################################
const std::string PROT_1 = "10100110100101100101";
const int FALTUNG_LENGTH = 20;
const int BEST_FITNESS = 999;

// Konfiguration
const int AMOUNT_PROTEINS = 10;
const int GENERATIONS = 100;
const int ELITE_AMOUNT = 3;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 50; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 50;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = true;
/*/

/*/### PROTEIN 1 ######################################################
const std::string PROT_1 = "10100110100101100101";
const int FALTUNG_LENGTH = 20;
const int BEST_FITNESS = 9;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 1000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 50;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/

/*### PROTEIN 2 ######################################################
const std::string PROT_1 = "110010010010010010010011";
const int FALTUNG_LENGTH = 24;
//const int BEST_FITNESS = 9;
const int BEST_FITNESS = 99;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/

/*### PROTEIN 3 ######################################################
const std::string PROT_1 = "0010011000011000011000011";
const int FALTUNG_LENGTH = 25;
const int BEST_FITNESS = 8;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/

/*### PROTEIN 4 ######################################################
const std::string PROT_1 = "000110011000001111111001100001100100";
const int FALTUNG_LENGTH = 36;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 13; // Meistens
//const int BEST_FITNESS = 14;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 200; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 8; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/

/*### PROTEIN 5 ######################################################
const std::string PROT_1 = "001001100110000011111111110000001100110010011111";
const int FALTUNG_LENGTH = 48;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 22;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 0;
const double MUTATION_RATE_VARIANCE = 4; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/

/*### PROTEIN 6 ######################################################
const std::string PROT_1 = "11010101011110100010001000010001000101111010101011";
const int FALTUNG_LENGTH = 50;
const int BEST_FITNESS = 99;
//const int BEST_FITNESS = 21;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 10; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/
//*### PROTEIN 100er ######################################################
const std::string PROT_1 = "0001100111100111011011011110000000011111100111111000000000101101111111111100111011010010111000000111";
const int FALTUNG_LENGTH = 100;
const int BEST_FITNESS = 50;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = 100;
const double GENERATION_WITH_LOWEST_MUTATION_RATE = 1000; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 40;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 10; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/
#pragma endregion

// Elite Auswahl wird nur einmal übernommen, danach nur dessen Mutationen
const bool ELITE_ONLY_ONCE = false;

// Doppelete erlauben
const bool ALLOW_DOUBLES = false;

// Die Senkung der Mutationsrate pro Generation
const double MUTATION_RATE_LOWER = (MUTATION_RATE - MUTATION_RATE_LOWEST) / GENERATION_WITH_LOWEST_MUTATION_RATE;

// Anzeige
int NUMBER_WIDTH = 3;

// Für Farbausgabe
HANDLE hstdout = GetStdHandle(STD_OUTPUT_HANDLE);

void main()
{
	#pragma region Initialisierung und mehr
	time_t t;
    time(&t);
    srand((unsigned int)t);

	std::list<Protein*> populations = std::list<Protein*>();
	
	Protein bestProtein = Protein(0);
	int bestProtein_round = -1;
	Protein* highscore[AMOUNT_PROTEINS];
	std::list<Protein*>::iterator populationIter;
	double generationFitnessTotal;
	double generationFitnessAverage;
	int proteinIndex;

	// Ausgangsprotein
	Protein* proteinToCopy = new Protein(PROT_1);

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
	if(SHOW_DETAILS) std::cout << std::setw(6) << std::right << "Gen    1 -> ";
	for(int i = 0; i<AMOUNT_PROTEINS; i++)
	{
		// Ausgangsprotein kopieren
		Protein* addMe = new Protein(FALTUNG_LENGTH);
		addMe = proteinToCopy->Copy();

		if(SHOW_DETAILS)
		{
			int f = addMe->CalcFitness();
			
			std::cout << std::setw(NUMBER_WIDTH) << std::right << f;
		}

		populations.push_front(addMe);
	}

	std::cout << std::endl;
	#pragma endregion
	
	#pragma region In allen Generationen, alle Proteine durchgehen
	for(int i = 0; i<GENERATIONS; i++)
	{
		bool showRow = (SHOW_DETAILS || i%50 == 0);

		if(showRow) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> ";

		#pragma region Mutieren
		// Mutieren lassen und Besten suchen
		generationFitnessTotal = 0;
		proteinIndex = 0;
		for (populationIter = populations.begin(); populationIter != populations.end(); ++populationIter)
		{
			double mutation_rate = MUTATION_RATE;

			// Toleranz in Mutationsrate einbauen
			double tolerance = 0;
			if(MUTATION_RATE_VARIANCE != 0)
			{
				// Auf Mutation die Toleranz dazurechnen
				tolerance = (rand() % (int)(MUTATION_RATE_VARIANCE * 10000) - (MUTATION_RATE_VARIANCE/2.0 * 10000)) / 10000.0;

				mutation_rate += tolerance;
				
				// Auf Grenzen achten
				if(mutation_rate < MUTATION_RATE_LOWEST) mutation_rate = MUTATION_RATE_LOWEST;
			}

			// Mutieren (aber nicht die "Elite", die bleibt bestehen)
			if(proteinIndex >= ELITE_AMOUNT) (*populationIter)->Mutate(mutation_rate);

			// Fitness ausrechnen
			int f = (*populationIter)->CalcFitness();

			// Fitness in die Durschnittsrechnung einrechnen
			generationFitnessTotal += (f == -1 ? 0 : f);

			// Die Elite nicht in den Highscore (und damit nicht in die nächste Generation übernehemen)
			if(!ELITE_ONLY_ONCE || proteinIndex >= ELITE_AMOUNT) highscore[proteinIndex - (ELITE_ONLY_ONCE ? ELITE_AMOUNT : 0)] = *populationIter;

			// Fortlaufender index;
			proteinIndex++;
		}

		// Highscore sortieren
		//sortHighscore(AMOUNT_PROTEINS - (ELITE_ONLY_ONCE ? ELITE_AMOUNT : 0), highscore);
		qsort (highscore, AMOUNT_PROTEINS - (ELITE_ONLY_ONCE ? ELITE_AMOUNT : 0), sizeof(Protein*), compare);

		// Bestes Protein merken
		Protein bestGenerationProtein = *highscore[0]->Copy();

		// Durschnittsrechnung abschließen
		generationFitnessAverage = ceil(generationFitnessTotal / populations.size() * 10) / 10.0;
		#pragma endregion

		#pragma region Anzeigen
		if(SHOW_DETAILS)
		{
			for (populationIter = populations.begin(); populationIter != populations.end(); ++populationIter)
			{
				int f = (*populationIter)->CalcFitness();

				if(f == -1) SetConsoleTextAttribute(hstdout, 0x0c);
				else if((*populationIter)->Identifier() == bestGenerationProtein.Identifier()) SetConsoleTextAttribute(hstdout, 0xa0);

				if(showRow) std::cout << std::setw(NUMBER_WIDTH) << std::right << f;
			
				SetConsoleTextAttribute(hstdout, 0x0f);
			}
		}
		#pragma endregion
		
		#pragma region Auswerten, welches bestes Protein ist
		bool isNewBest = (bestGenerationProtein.CalcFitness() > bestProtein.CalcFitness()); 
		
		if(isNewBest)
		{
			bestProtein = bestGenerationProtein;
			
			bestProtein_round = i + 2;
		}

		// Auswertung anzeigen
		if(!showRow && isNewBest) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> "; // Wenn SUCCESSES_ONLY, dann hier einen Ersatz anzeigen

		if(showRow || isNewBest)
		{
			std::cout << " |" << std::setw(5) << generationFitnessAverage;
		
			std::cout << std::setw(5) << bestGenerationProtein.CalcFitness();
		
			if(isNewBest) SetConsoleTextAttribute(hstdout, 0xa0);
			std::cout << std::setw(5) << bestProtein.CalcFitness();
			SetConsoleTextAttribute(hstdout, 0x0f);

			std::cout << std::setw(10) << MUTATION_RATE;

			std::cout << std::endl;
		}
		#pragma endregion
		
		#pragma region Selektion
		//*
		// X Elite-Proteine wählen, die nicht gleich sind
		std::list<std::string> proteinIds = std::list<std::string>();
		std::list<Protein*> newPopulation = std::list<Protein*>();
		for(int i = 0; i<AMOUNT_PROTEINS - ELITE_AMOUNT; i++)
		{
			// Id des Proteins bilden
			std::string id = highscore[i]->Identifier();
			
			// Wenn Protein noch nicht in Liste vorhanden, dann hinzufügen
			if(!ALLOW_DOUBLES)
			{
				std::list<std::string>::iterator findIter = std::find(proteinIds.begin(), proteinIds.end(), id);
				if(findIter == proteinIds.end())
				{
					newPopulation.push_back(highscore[i]->Copy());
					proteinIds.push_back(id);
				}

				// Wenn Liste voll, dann abbrechen
				if(proteinIds.size() == ELITE_AMOUNT) break;
			}
			else newPopulation.push_back(highscore[i]->Copy());
		}

		// Wenn zu wenig Elite Proteine vorhanden, dann Rest mit Bestem Protein auffüllen
		if(proteinIds.size() < ELITE_AMOUNT)
		{
			for(int i = proteinIds.size(); i<ELITE_AMOUNT; i++) newPopulation.push_back(highscore[0]->Copy());
		}

		// Wenn Maximum erreicht, dann abbrechen
		if(bestProtein.CalcFitness() == BEST_FITNESS) break;

		// Neue Population generieren
		int proteinsLeft = AMOUNT_PROTEINS - ELITE_AMOUNT;
		for (populationIter = populations.begin(); populationIter != populations.end(); ++populationIter)
		{
			int f = (*populationIter)->CalcFitness();
			if(f == -1) f = 0;

			int amountOfClones = floor((f / generationFitnessTotal) * (AMOUNT_PROTEINS - ELITE_AMOUNT)); // floor(x + 0.5) simuliert round

			for(int a = 0; a < amountOfClones && proteinsLeft > 0; a++)
			{
				newPopulation.push_back((*populationIter)->Copy());
			
				proteinsLeft--;
			}

			if(proteinsLeft == 0) break;
		}

		// Alte Population löschen
		for (std::list<Protein*>::iterator deleteIter = populations.begin(); deleteIter != populations.end(); ++deleteIter) delete *deleteIter;

		populations = newPopulation;
		#pragma endregion

		// Mutationsrate senken
		MUTATION_RATE -= MUTATION_RATE_LOWER;

		// Wenn niedrigste Mutationsrate erreicht, dann Untergrenze setzen
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
		int f = highscore[i]->CalcFitness();
		std::string id = highscore[i]->Identifier();
		std::cout << "Rank " << (i+1) << " (Fitness: " << f << ", Id: " << id << ")" << std::endl;
		highscore[i]->Show();
	}
	#pragma endregion
}