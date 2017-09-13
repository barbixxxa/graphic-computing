#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
	//Definir altura e largura da imagem
	float altura = 600;
	float largura = 800;
	
	ofstream outFile;
	//Criar o arquivo .pgm
	outFile.open("arquivo.pgm");

	//Cabeçalho do arquivo
	outFile << "P3" << endl; //P3 = RGB e P2 = Cinza
	outFile << "800 600" << endl; // Resulucao da imagem
	outFile << "255" << endl; // Intensidade maxima do pixel

	//Varrer a imagem linha por linha e coluna por coluna
	for (int linha = 0; linha < altura; ++linha) {
		for (int coluna = 0; coluna < largura; ++coluna) {
			if (coluna < (largura / 3)) { //Pixel azul
				outFile << "0 85 164 ";
				
			}else if(coluna > (largura/3) && coluna < (largura * 2) / 3){ //Pixel Branco
				outFile << "255 255 255 "; 
			}
			else { //Pixel vermelho
				outFile << "239 65 53 ";

			}
		}
	}
	
	
	outFile.close(); //Fechar o arquivo
	return 0;
}
