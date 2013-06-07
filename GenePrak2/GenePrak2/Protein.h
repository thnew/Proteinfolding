#pragma once
#include <iostream>

class Amino;

class Protein
{
	private:
		Amino* Aminos;
		int Length;
		int FamilyNr;
		bool NoMutation;

		// Caching Attribute
		int Fitness;
		std::string Id;

	public:
		Protein(int length, int familyNr);
		Protein(std::string proteinString, int familyNr);
		Protein();
		~Protein();

		Amino* GetAminos();
		void SetNoMutation(bool);

		int CalcFitness();
		Protein* Copy();
		std::string Identifier();
		void Mutate(double mutationRate);
		void Show();
};