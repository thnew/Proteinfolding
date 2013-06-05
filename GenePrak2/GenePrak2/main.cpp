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

void sortHighscore(int, Protein**);

#pragma region Proteine
/*### Testprotein ######################################################
const std::string PROT_1 = "10100110100101100101";
const int FALTUNG_LENGTH = 20;
const int BEST_FITNESS = 999;

// Konfiguration
const int AMOUNT_PROTEINS = 10;
const int GENERATIONS = 100;
const int ELITE_AMOUNT = 3;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 50; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 50;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 5; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = true;
//*/

//*### PROTEIN 1 ######################################################
const std::string PROT_1 = "10100110100101100101";
const int FALTUNG_LENGTH = 20;
const int BEST_FITNESS = 9;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 1000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
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
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
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
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
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
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 200; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
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
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
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
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 100; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 75;
const double MUTATION_RATE_LOWEST = 10;
const double MUTATION_RATE_VARIANCE = 10; // Toleranz der Mutationsrate/Um wieviel die Rate beim Mutieren maximal abweichen darf

const bool SHOW_DETAILS = false;
//*/
/*### PROTEIN 100er ######################################################
const std::string PROT_1 = "0001100111100111011011011110000000011111100111111000000000101101111111111100111011010010111000000111";
const int FALTUNG_LENGTH = 100;
const int BEST_FITNESS = 50;

// Konfiguration
const int AMOUNT_PROTEINS = 500;
const int GENERATIONS = 5000;
const int ELITE_AMOUNT = AMOUNT_PROTEINS / 20;
const double GENERATION_WITH_LOWEST_MUATTION_RATE = 500; // Die Generation, bei der die niedrigste Mutationsrate erreicht werden soll
double MUTATION_RATE = 50;
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
const double MUTATION_RATE_LOWER = (MUTATION_RATE - MUTATION_RATE_LOWEST) / GENERATION_WITH_LOWEST_MUATTION_RATE;

// Anzeige
int NUMBER_WIDTH = 3;

// Für Farbausgabe
HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );

void main()
{
	#pragma region Initialisierung und mehr
	time_t t;
    time(&t);
    srand((unsigned int)t);

	std::list<Protein*> populations = std::list<Protein*>();
	Protein* eliteProteins[ELITE_AMOUNT];

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
	Protein* bestProtein = new Protein(0);
	int bestProtein_round= -1;
	Protein* highscore[AMOUNT_PROTEINS];
	for(int i = 0; i<GENERATIONS; i++)
	{
		bool showRow = SHOW_DETAILS || i%500 == 0;

		if(showRow) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> ";

		#pragma region Mutieren

		// Mutieren lassen und Besten suchen
		std::list<Protein*>::iterator iter;
		double generationAverage = 0;
		int count = 0;
		for (iter = populations.begin(); iter != populations.end(); ++iter)
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
			if(count >= ELITE_AMOUNT) (*iter)->Mutate(mutation_rate);

			// Fitness ausrechnen
			int f = (*iter)->CalcFitness();

			// Fitness in die Durschnittsrechnung einrechnen
			generationAverage += (f == -1 ? 0 : f);

			// Die Elite nicht in den Highscore (und damit in die nächste Generation übernehemen)
			if(!ELITE_ONLY_ONCE || count >= ELITE_AMOUNT) highscore[count - (ELITE_ONLY_ONCE ? ELITE_AMOUNT : 0)] = *iter;

			// Fortlaufender index;
			count++;
		}

		// Highscore sortieren
		sortHighscore(AMOUNT_PROTEINS - (ELITE_ONLY_ONCE ? ELITE_AMOUNT : 0), highscore);

		// Auswerten, welches bestes Protein ist
		Protein* bestGenerationProtein = highscore[0];

		// Durschnittsrechnung abschließen
		generationAverage = ceil(generationAverage / populations.size() * 10) / 10.0;
		#pragma endregion

		#pragma region Anzeigen
		if(SHOW_DETAILS)
		{
			for (iter = populations.begin(); iter != populations.end(); ++iter)
			{
				int f = (*iter)->CalcFitness();

				if(f == -1) SetConsoleTextAttribute(hstdout, 0x0c);
				else if(*iter == bestGenerationProtein) SetConsoleTextAttribute(hstdout, 0xa0);

				if(showRow) std::cout << std::setw(NUMBER_WIDTH) << std::right << f;
			
				SetConsoleTextAttribute(hstdout, 0x0f);
			}
		}
		#pragma endregion
		
		#pragma region Auswerten, welches bestes Protein ist
		bool isNewBest = (bestGenerationProtein->CalcFitness() > bestProtein->CalcFitness()); 
		
		if(isNewBest)
		{
			bestProtein = bestGenerationProtein;
			
			bestProtein_round = i + 2;
		}

		// Auswertung anzeigen
		if(!showRow && isNewBest) std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> "; // Wenn SUCCESSES_ONLY, dann hier einen Ersatz anzeigen

		if(showRow || isNewBest)
		{
			std::cout << " |" << std::setw(5) << generationAverage;
		
			std::cout << std::setw(5) << bestGenerationProtein->CalcFitness();
		
			if(isNewBest) SetConsoleTextAttribute(hstdout, 0xa0);
			std::cout << std::setw(5) << bestProtein->CalcFitness();
			SetConsoleTextAttribute(hstdout, 0x0f);

			std::cout << std::setw(10) << MUTATION_RATE;

			std::cout << std::endl;
		}
		#pragma endregion
		
		#pragma region Selektion
		// X Elite-Proteine wählen, die nicht gleich sind
		std::list<std::string> proteinIds = std::list<std::string>();
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
					eliteProteins[proteinIds.size()] = highscore[i];
					proteinIds.push_back(id);
				}

				// Wenn Liste voll, dann abbrechen
				if(proteinIds.size() == ELITE_AMOUNT) break;
			}
			else eliteProteins[proteinIds.size()] = highscore[i];
		}

		// Wenn zu wenig Elite Proteine vorhanden, dann Rest mit Bestem Protein auffüllen
		if(proteinIds.size() < ELITE_AMOUNT)
		{
			for(int i = proteinIds.size(); i<ELITE_AMOUNT; i++) eliteProteins[i] = highscore[0];
		}

		// Neue Population generieren
		std::list<Protein*> newPopulation = std::list<Protein*>();
		for(int i = 0; i<AMOUNT_PROTEINS; i++)
		{
			// Wenn es eine Elite gibt, dann aus dieser die zu Übernehmenden wählen
			int copyFrom = (ELITE_AMOUNT != 0 ? (i % ELITE_AMOUNT) : i);

			// Die Elite Proteine werden unverändert weitergegeben
			if(i < ELITE_AMOUNT) newPopulation.push_back(eliteProteins[copyFrom]);
			// Die restlichen werden kopiert und gespeichert
			else newPopulation.push_back(eliteProteins[copyFrom]->Copy());
		}

		populations = newPopulation;
		#pragma endregion

		// Wenn Maximum erreicht, dann abbrechen
		if(bestProtein->CalcFitness() == BEST_FITNESS) break;

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
		int f = eliteProteins[i]->CalcFitness();
		std::string id = eliteProteins[i]->Identifier();
		std::cout << "Rank " << (i+1) << " (Fitness: " << f << ", Id: " << id << ")" << std::endl;
		eliteProteins[i]->Show();
	}
	#pragma endregion
}

void sortHighscore(int length, Protein** sortMe)
{
	Protein* temp_protein;
   
	for(int i=0; i<length; i++)
	{
		for(int j=0; j<length-1; j++)
		{
			if(sortMe[j]->CalcFitness() < sortMe[j+1]->CalcFitness())
			{
				temp_protein = sortMe[j];
				sortMe[j] = sortMe[j+1];
				sortMe[j+1] = temp_protein;
			}
		}
	}
}