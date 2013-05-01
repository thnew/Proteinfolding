#pragma once

enum SkyDirection {NORTH = 0, EAST = 1, SOUTH= 2, WEST= 3};
enum Direction { FORWARD = 888, LEFT = 988, RIGHT = 889 };

class Amino
{
	public:
		class DisplayInfo
		{
			private:
				bool Set;
				SkyDirection SkyDir;
				Direction RelDir;
				bool Hydrophob;
			public:
				DisplayInfo(){ this->Set = false; }
				DisplayInfo(SkyDirection sDir, bool hydrophob)
				{
					this->SkyDir = sDir;
					this->Set = true;
					this->Hydrophob = hydrophob;
				}
				~DisplayInfo(){}
				
				void SetSkyDir(SkyDirection sDir)
				{
					this->SkyDir = sDir;
				}void SetRelDir(Direction dir)
				{
					this->RelDir = dir;
				}
				void SetHydrophob(bool isHydrophob)
				{
					this->Hydrophob = isHydrophob;
				}
				
				SkyDirection GetSkyDir()
				{
					return this->SkyDir;
				}
				Direction GetRelDir()
				{
					return this->RelDir;
				}
				bool IsHydrophob()
				{
					return this->Hydrophob;
				}
				bool IsSet()
				{
					return this->Set;
				}
		};
	private:
		Direction Dir;
		bool Hydrophob;
	public:
		Amino(Direction, bool);
		Amino(void);
		~Amino(void);

		// Getter
		Direction GetDir();
		bool IsHydrophob();

		// Setter
		void SetDir(Direction);
		void SetHydrophob(bool);
};

