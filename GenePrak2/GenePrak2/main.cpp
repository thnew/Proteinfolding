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
void mutate(Amino* faltung);

/*
int GetMinX(int*);
int GetMaxX(int*);
int GetMinY(int*);
int GetMaxY(int*);
//*/

const int FALTUNG_LENGTH = 30;

// Für Farbausgabe
HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );

void main()
{
	Amino* faltung = new Amino[FALTUNG_LENGTH];
	
	int amountProteins = 15;
	int generations = 200;

	std::list<Amino*> populations = std::list<Amino*>();

	int numberWidth = 3;

	// Titel
	SetConsoleTextAttribute(hstdout, 0x1f);
	std::cout << "Proteins -> ";
	for(int i = 0; i<amountProteins; i++) std::cout << std::setw(numberWidth) << std::right << i;
	SetConsoleTextAttribute(hstdout, 0x0f);
	std::cout << std::endl;

	// Population generieren
	std::cout << std::setw(6) << std::right << "Gen    1 -> ";
	for(int i = 0; i<amountProteins; i++)
	{
		Amino* addMe = new Amino[FALTUNG_LENGTH];

		int f = fitness(addMe);

		std::cout << std::setw(numberWidth) << std::right << f;

		populations.push_front(addMe);
	}

	std::cout << std::endl;

	// Alle Proteine durchgehen
	Amino* bestProtein = new Amino[FALTUNG_LENGTH];
	int bestProtein_fitness = -1;
	for(int i = 0; i<generations; i++)
	{
		std::cout << "Gen "<< std::setw(4) << std::right << (i+2) << " -> ";

		// Mutieren lassen und Besten suchen
		std::list<Amino*>::iterator iter;
		Amino* bestGenerationProtein = 0;
		int bestGenerationProtein_fitness = -1;
		for (iter = populations.begin(); iter != populations.end(); ++iter)
		{
			mutate(*iter);

			int f = fitness(*iter);

			if(f > bestGenerationProtein_fitness)
			{
				bestGenerationProtein_fitness = f;
				bestGenerationProtein = *iter;
			}
		}

		// Anzeigen
		for (iter = populations.begin(); iter != populations.end(); ++iter)
		{
			int f = fitness(*iter);

			if(f == -1) SetConsoleTextAttribute(hstdout, 0x0c);
			else if(*iter == bestGenerationProtein) SetConsoleTextAttribute(hstdout, 0xa0);

			std::cout << std::setw(numberWidth) << std::right << f;
			
			SetConsoleTextAttribute(hstdout, 0x0f);
		}

		if(bestGenerationProtein_fitness > bestProtein_fitness)
		{
			//bestProtein = bestGenerationProtein;
			memcpy(bestProtein, bestGenerationProtein, FALTUNG_LENGTH * sizeof(*bestGenerationProtein));
			bestProtein_fitness = bestGenerationProtein_fitness;
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;

	std::cout << "BEST PROTEIN with fitness " << bestProtein_fitness << std::endl;
	std::cout << "Fitness recalculated: " << fitness(bestProtein) << std::endl;

	show(bestProtein);

	/*
	faltung[0] = Amino(FORWARD, false);
	faltung[1] = Amino(FORWARD, true);
	faltung[2] = Amino(RIGHT, false);
	if(FALTUNG_LENGTH >= 4) faltung[3] = Amino(RIGHT, true);
	if(FALTUNG_LENGTH >= 5) faltung[4] = Amino(FORWARD, false);
	if(FALTUNG_LENGTH >= 6) faltung[5] = Amino(FORWARD, true);
	if(FALTUNG_LENGTH >= 7) faltung[6] = Amino(RIGHT, false);
	//*/

	/*
	Amino* lastFitness = new Amino[FALTUNG_LENGTH];
	memcpy(lastFitness, faltung, sizeof(faltung));

	int bestFitness = -1;

	for(int i = 0; i<500; i++)
	{
		mutate(faltung);
		
		int f = fitness(faltung);
		
		if(f == -1) std::cout << "Invalid";
		else if(f == -2) std::cout << "Algorithm calculated a fitness smaller than 0, what is not possible";
		else std::cout << "Fitness: " << f;

		// Wenn Fitness schlechter als die der besten Faltung, dann die Beste laden
		if(f < bestFitness)
		{
			memcpy(faltung, lastFitness, sizeof(lastFitness));
			
			std::cout << " -> not good enough";
		}
		else
		{
			memcpy(lastFitness, faltung, sizeof(faltung));
			bestFitness = f;
		}

		std::cout << std::endl;
	}

	std::cout << std::endl << "BestFitness: " << bestFitness;
	//*/

	std::cout << std::endl;

	//system("PAUSE");
}

void mutate(Amino* faltung)
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

int fitness(Amino* faltung)
{
	int startX = 100, startY = 100;
	int x, y;
	int width = 200;
	int height = 200;
	int matrix[200][200];

	int connections = 0;
	int directConnections = 0;
	bool lastIsHydrophob = false;

	bool show = false;

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
		//*/
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

			if(show) std::cout << (amino->GetDir() == FORWARD ? "FORWARD" : (amino->GetDir() == LEFT ? "LEFT" : "RIGHT")) << " | ";
			
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

			/*/ zusätzlich eintragen, ob Nord oder Süd
			matrix[x][y] |= (amino->GetDir() == NORTH) * 10000;
			matrix[x][y] |= (amino->GetDir() == EAST) * 1000;
			matrix[x][y] |= (amino->GetDir() == SOUTH) * 100;
			matrix[x][y] |= (amino->GetDir() == WEST) * 10;
			//*/

			//std::cout << x << "/" << y << ": " << amino->IsHydrophob() << std::endl;
			
			lastIsHydrophob = amino->IsHydrophob();
		}
		#pragma endregion

		#pragma region Matrix ausgeben
		if(show)
		{
			/*/ Die Maße des Protins errechnen
			int minX = GetMinX(matrix);
			int maxX = GetMaxX(matrix);
			int minY = GetMinY(matrix);
			int maxY = GetMaxY(matrix);
			//*/
			
			const int protWidth = maxX - minX;
			const int protHeight = maxY - minY;

			std::cout << "Protein: (X: " << minX << "/" << maxX << " | Y: " << minY << "/" << maxY << ")" << std::endl;

			for(int y = maxY; y>=minY; y--)
			{
				for(int x = minX; x<=maxX; x++)
				{
					int hydrophob = matrix[x][y];

					// Wenn leer, dann weiter
					if(hydrophob == -1) std::cout << ' ';
					else if(hydrophob == 0) std::cout << '0';
					else if(hydrophob == 1) std::cout << 'X';
				}

				std::cout << std::endl;
			}
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
					//std::cout << x << "/" << y << std::endl;
					
					//std::cout << "Before:" << connections;
					
					// Alle Verbindungen des Hydrophobs zählen
					if(matrix[x+1][y] == 1) connections++;
					if(matrix[x-1][y] == 1) connections++;
					if(matrix[x][y+1] == 1) connections++;
					if(matrix[x][y-1] == 1) connections++;

					//std::cout << "After:" << connections << std::endl;

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
	int startX = 100, startY = 100;
	int x, y;
	int width = 200;
	int height = 200;
	Amino::DisplayInfo* matrix[200][200];

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
	/*/ Die Maße des Protins errechnen
	int minX = GetMinX(matrix);
	int maxX = GetMaxX(matrix);
	int minY = GetMinY(matrix);
	int maxY = GetMaxY(matrix);
	//*/
			
	//const int protWidth = maxX - minX;
	//const int protHeight = maxY - minY;

	maxX++;
	maxY++;
	minX--;
	minY--;

	std::cout << "Protein: (X: " << minX << "/" << maxX << " | Y: " << minY << "/" << maxY << ")" << std::endl;

	for(int y = maxY; y>=minY; y--)
	{
		for(int x = minX; x<=maxX; x++)
		{
			Amino::DisplayInfo* amino = matrix[x][y];

			// Wenn leer, dann weiter
			if(amino == nullptr) std::cout << ' ';
			else
			{
				if(!amino->IsHydrophob()) SetConsoleTextAttribute(hstdout, 0xcf);
				else if(amino->IsHydrophob()) SetConsoleTextAttribute(hstdout, 0xaf);

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

/*
int GetMinX(Amino* faltung)
{
	int minX = 0;
	int x = 0;

	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		if(faltung[i].GetDir() == EAST) x++;
		else if(faltung[i].GetDir() == WEST) x--;

		if(x < minX) minX = x;
	}

	return minX;
}

int GetMaxX(Amino* faltung)
{
	int maxX = 0;
	int x = 0;

	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		if(faltung[i].GetDir() == EAST) x++;
		else if(faltung[i].GetDir() == WEST) x--;

		if(x > maxX) maxX = x;
	}

	return maxX;
}

int GetMinY(Amino* faltung)
{
	int minY = 0;
	int y = 0;

	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		if(faltung[i].GetDir() == NORTH) y++;
		else if(faltung[i].GetDir() == SOUTH) y--;

		if(y < minY) minY = y;
	}

	return minY;
}

int GetMaxY(Amino* faltung)
{
	int maxY = 0;
	int y = 0;

	for(int i = 0; i<FALTUNG_LENGTH; i++)
	{
		if(faltung[i].GetDir() == NORTH) y++;
		else if(faltung[i].GetDir() == SOUTH) y--;

		if(y > maxY) maxY = y;
	}

	return maxY;
}
//*/