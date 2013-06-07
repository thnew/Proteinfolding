#include "Protein.h"
#include "Amino.h"
#include <windows.h>
#include <iomanip>
#include <windows.h>

Protein::Protein(int length, int familyNr)
{
	this->Fitness = -5;
	this->Length = length;
	this->Aminos = new Amino[length];
	this->FamilyNr = familyNr;
	this->NoMutation = false;
}

Protein::Protein(std::string copyMe, int familyNr)
{
	this->Fitness = -5;
	this->Length = copyMe.size();
	this->Aminos = new Amino[this->Length];
	this->FamilyNr = familyNr;
	this->NoMutation = false;

	for(int i = 0; i<this->Length; i++) this->Aminos[i] = Amino(FORWARD, copyMe[i] == 49);
}

Protein::Protein()
{
	this->Fitness = -5;
	this->NoMutation = false;
}

Protein::~Protein()
{
	//delete this->Aminos;
}

Amino* Protein::GetAminos()
{
	return this->Aminos;
}

int Protein::CalcFitness()
{
	if(this->Fitness != -5) return this->Fitness;

	//if(this->Length < 0) std::cout << std::endl << "Fitness-> Length < 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	
	const int MATRIX_SIZE = this->Length * 2 + 3;
	int x, y;
	const int width = MATRIX_SIZE;
	const int height = MATRIX_SIZE;
	int** matrix = new int*[width];
	int startX = (width - 1) / 2;
	int startY = (height - 1) / 2;

	int connections = 0;
	int directConnections = 0;
	bool lastIsHydrophob = false;

	try
	{
		#pragma region Initialisierung
		// Alle Felder mit -1 (==leer) initialisieren
		for(int pX = 0; pX<height; pX++)
		{
			matrix[pX] = new int[width];

			for(int pY = 0; pY<width; pY++)
			{
				matrix[pX][pY] = -1;
			}
		}

		// Anfangskoordinaten setzen
		x = startX;
		y = startY;
		#pragma endregion

		#pragma region Aminos in Matrix eintragen
		int dir = NORTH;
		//int maxX = startX, minX = startX, maxY = startY, minY= startY;
		for(int i = 0; i<this->Length; i++)
		{
			Amino* amino = &this->Aminos[i];

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

			/*/ Grenzen ermitteln
			if(x > maxX) maxX = x;
			if(x < minX) minX = x;
			if(y > maxY) maxY = y;
			if(y < minY) minY = y;
			//*/

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
		// Matrix löschen
		for(int pX = 0; pX<height; pX++) delete matrix[pX];
		delete matrix;

		this->Fitness = e;
		return e;
	}

	// Matrix löschen
	for(int pX = 0; pX<height; pX++) delete matrix[pX];
	delete matrix;

	// Alle drekten Verbundungen von der Gesamtzahl aller Verbindungen abziehen
	int f = connections - directConnections;

	// Wenn die gezählte Zahl kleiner 0 ist, dann liegt ein Fehler vor, der mit einer -2 signalisiert wird
	if(f < 0) return -2;

	// Fitness merken
	this->Fitness = f;

	// Fitness zurückgeben
	return f;
}

void Protein::SetNoMutation(bool noMutation)
{
	this->NoMutation = noMutation;
}

Protein* Protein::Copy()
{
	Protein* returnMe = new Protein(this->Length, this->FamilyNr);

	returnMe->SetNoMutation(this->NoMutation);

	memcpy(returnMe->GetAminos(), this->Aminos, this->Length * sizeof(*this->Aminos));

	return returnMe;
}

void Protein::Mutate(double mutationRate)
{
	if(this->NoMutation)
	{
		this->NoMutation = false;
		return;
	}

	//if(this->Length < 0) std::cout << std::endl << "Mutate-> Length < 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	//std::cout << "Rate: " << mutationRate << std::endl;
	
	// Anzahl der Mutationen berechnen
	int mutations = (mutationRate / 100.0) * this->Length;

	// Verhindern, dass Mutationsrate so klein ist, dass gar keine Mutataion mehr stattfindet
	mutations = (mutations < 1 ? 1 : mutations);

	for(int i=0; i<mutations; i++)
	{
		// Nach Zufall eine Aminosäure ermitteln
		int randAmino = rand() % (this->Length - 1) + 1; // Erste Aminosäure ignorieren
	
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
		} while(randDir == this->Aminos[randAmino].GetDir());

		// Richtung in der Aminosäure setzen
		this->Aminos[randAmino].SetDir((Direction)randDir);
	}

	// Fitness zurücksetzen
	this->Fitness = -5;

	// Id zurücksetzen
	this->Id = "";
}

// Generates a number, that identifies a protein
const char hexchar[17] = "0123456789ABCDEF";
std::string Protein::Identifier()
{
	/*
	int id = 0;

	for(int i = 0; i<this->Length; i++)
	{
		int aminoId = 0;
		
		if(this->Aminos[i].GetDir() == LEFT) aminoId = 0;
		else if(this->Aminos[i].GetDir() == FORWARD) aminoId = 1;
		else if(this->Aminos[i].GetDir() == RIGHT) aminoId = 2;
		
		///id += aminoId * (i+1) * 100;

		id = id << 2;
		id += aminoId;
	}

	return id;
	/*/

	if(this->Id.size() > 0) return this->Id;
	
	int familyIdLength = 4;

	//std::string id = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
	//std::string id = "000000000000000000000000000000";
	std::string id = new char[this->Length];
	id = id.substr(0,this->Length + familyIdLength);

	// FamilyId ermitteln und einfügen
	int familyNr = this->FamilyNr;
	for(int i = 0; i < familyIdLength; i++)
	{
		id[i] = (familyNr > 0 ? hexchar[familyNr & 0xf] : '0');
		familyNr = familyNr >> 4;
	}

	id[familyIdLength] = '-';

	// restliche Id bilden
	char aminoId = '0';
	Direction dir;
	for(int i = familyIdLength + 1; i<this->Length + familyIdLength; i++)
	{
		dir = this->Aminos[i - familyIdLength].GetDir();

		if(dir == LEFT) aminoId = 'L';
		else if(dir == FORWARD) aminoId = 'F';
		else if(dir == RIGHT) aminoId = 'R';
		
		id[i] = aminoId;
	}

	this->Id = id;

	return id;
	//*/
}

void Protein::Show()
{
	// Für Farbausgabe
	HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );
	
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
	for(int i = 0; i<this->Length; i++)
	{
		Amino* amino = &this->Aminos[i];

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