#include <fstream>
#include <armadillo>
#include <thread>
#include <cmath>

using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html


bool intersect(const vec& origin, const vec& direction, const vec& center, double radius, vec& normal, vec& hit) {
	vec oc = origin - center;
	double a = dot(direction, direction);
	double b = 2.0*dot(direction, oc);
	double c = dot(oc, oc) - radius*radius;

	double delta = b*b - 4.0*a*c;
	if (delta < 0.0) {
		return false;
	}
	double t1 = (-b + sqrt(delta)) / (2.0*a);
	double t2 = (-b - sqrt(delta)) / (2.0*a);
	double t = (t1 < t2) ? t1 : t2;

	hit = origin + direction*t;
	normal = hit - center;
	return true;
}

struct Material {
	vec kd;
	vec ke;

	double shininess;
	double reflection;

	Material& operator = (const Material& mat) {
		ke = mat.ke;
		kd = mat.kd;
		shininess = mat.shininess;
		reflection = mat.reflection;
		return *this;
	}
};

struct Ambient {
	vec ka;

	Ambient& operator = (const Ambient& amb) {
		ka = amb.ka;
		return *this;
	}
};

struct Luz {
	vec pos;
	vec radiancia;

	Luz& operator = (const Luz& luz) {
		pos = luz.pos;
		radiancia = luz.radiancia;
		return *this;
	}
};

struct Intersection {
	double t;
	vec p;
	vec normal;
	Material material;
	bool flipNormal;

	Intersection(double t, const vec& p, const vec& normal, Material material) : t(t), p(p), normal(normal), material(material), flipNormal(false) {}
	Intersection() : t(-1.0), flipNormal(false) {}
	Intersection(const Intersection& i) : t(i.t), p(i.p), normal(i.normal), material(i.material), flipNormal(i.flipNormal) {}

	Intersection& operator = (const Intersection& i) {
		t = i.t;
		p = i.p;
		normal = i.normal;
		material = i.material;
		flipNormal = i.flipNormal;
		return *this;
	}
};

struct Ray {
	vec origin;
	vec direction;
	Ray(const vec& origin, const vec& direction) : origin(origin), direction(direction) {}
};

class Object {
protected:
	Material material_;

public:
	Object(const Material(m)) : material_(m) {}
	virtual ~Object() {}
	virtual bool intersect(const Ray& r, Intersection& da) const = 0;
};

class Sphere : public Object {
private:
	double radius_;
	vec center_;

public:
	Sphere(const vec& center, double radius, const Material& materialSphere) : Object(materialSphere), center_(center), radius_(radius) {}

	virtual bool intersect(const Ray& r, Intersection& da) const {

		vec oc = r.origin - center_;
		double a = dot(r.direction, r.direction);
		double b = 2.0*dot(r.direction, oc);
		double c = dot(oc, oc) - radius_*radius_;

		double delta = b*b - 4.0*a*c;

		if (delta < 0.0) { return false; }

		double t1 = (-b + sqrt(delta)) / (2.0*a);
		double t2 = (-b - sqrt(delta)) / (2.0*a);

		da.t = (t1 < t2) ? t1 : t2;
		da.p = r.origin + r.direction *da.t;
		da.normal = da.p - center_;
		da.material = material_;


		return true;
	}
};

class Triangle : public Object {
private:
	vec A_, B_, C_;
public:
	Triangle(const vec& A, const vec& B, const vec& C, const Material& materialTriangle) : Object(materialTriangle), A_(A), B_(B), C_(C) {}

	virtual bool intersect(const Ray& r, Intersection& da) const {

		mat M, M2, M3, M4;

		//A
		M.insert_cols(0, A_ - B_);
		M.insert_cols(1, A_ - C_);
		M.insert_cols(2, r.direction);

		//beta
		M2.insert_cols(0, A_ - r.origin);
		M2.insert_cols(1, A_ - C_);
		M2.insert_cols(2, r.direction);

		//gamma
		M3.insert_cols(0, A_ - B_);
		M3.insert_cols(1, A_ - r.origin);
		M3.insert_cols(2, r.direction);


		double detM = det(M);
		double beta = det(M2) / (detM);
		double gamma = det(M3) / (detM);

		const vec& ba = B_ - A_;
		const vec& ca = C_ - A_;
		const vec& oa = r.origin - A_;
		da.normal = cross(ba, ca);
		da.t = (-(dot(da.normal, oa))) / (dot(r.direction, da.normal));


		if ((gamma<0) || (gamma>1)) return false;
		if ((beta<0) || (beta>(1.0 - gamma))) return false;

		da.p = r.origin + r.direction *da.t;
		da.material = material_;
		da.flipNormal = true;
		return true;
	}
};

struct Vertice {
	vec coordenada;

	Vertice& operator = (const Vertice& vertice) {
		coordenada = vertice.coordenada;
		return *this;
	}
};

struct Face {
	int a, b, c;

	Face& operator = (const Face& face) {
		a = face.a;
		b = face.b;
		c = face.c;
		return *this;
	}
};

int findNext(string s, int start, char c)
{
	string sub = s.substr(start, s.length());
	int r;
	r = sub.find(c);
	return r;
}

vector<Object*> loadObject(const string nomeArq, Material materialTriangle) { // recebe o nome do arquivo e o material do objeto

	int f1, f2, f3;
	vec v1, v2, v3;
	std::vector<Vertice> vertices;
	std::vector<Face> faces;
	Face face;
	Vertice v;
	std::vector<Object*> objetos;
	mat rotY, rotX;
	float mult_coord, somar_z;
	double sent = 0.707106; //seno de 45 graus
	double cost = 0.707106; //cos de 45 graus
	double senj = -0.5; //seno de 30 graus
	double cosj = 0.866025; //cos de 45 graus

	rotY << cost << 0.0 << -sent << endr //matriz de rotacao
		<< 0.0 << 1.0 << 0.0 << endr
		<< sent << 0.0 << cost;

	rotX << 1.0 << 0.0 << 0.0 << endr //matriz de rotacao
		<< 0.0 << cosj << -senj << endr
		<< 0.0 << senj << cosj;


	std::fstream myfile(nomeArq, std::fstream::in);
	if (myfile.is_open()) {
		std::string line;
		while (getline(myfile, line)) {
			std::stringstream stream(line, std::stringstream::in);
			std::string type;
			stream >> type;

			if (type == "#" || type == "") {
				continue;
			}
			else
				if (type == "m") {
					stream >> mult_coord; //valor para multiplicar todas as coordenadas
				}
				else
					if (type == "z") {
						stream >> somar_z; // valor para somar na coordenada z
					}
					else
						if (type == "v") {
							float x, y, z;
							stream >> x >> y >> z; //pega as coordenadas x,y,z do arquivo
							v.coordenada << x << y << z; //cria o vertice com as coordenadas

							v.coordenada = rotY * v.coordenada; //multiplica as coordenadas pela matriz de rotacao em torn ode Y
							v.coordenada = rotX * v.coordenada; //multiplica as coordenadas pela matriz de rotacao em torn ode X
							v.coordenada.at(2) = v.coordenada.at(2) + somar_z; // soma um valor a coordenada z para que possa ser vista
							v.coordenada *= mult_coord; //multiplica todas as coordenadas por um valor, para aumentar ou diminuir o tamanho da img
							vertices.push_back(v); //adiciona o vertice ao vector

						}
						else
							if (type == "vn") {
								float t, r, s;
								stream >> t >> r >> s;
								continue;

							}
							else
								if (type == "f") {
									int f, g, h;
									stream >> f >> g >> h;
									face.a = (f - 1);
									face.b = (g - 1);//cria a face com os vertices
									face.c = (h - 1);
									faces.push_back(face);
								}
								else {
									//não faz nada
								}


		}



		for (int i = 0; i < faces.size(); ++i) {

			f1 = faces[i].a; //pega os pontos que fazem parte da face
			f2 = faces[i].b;
			f3 = faces[i].c;

			v1 = vertices[f3].coordenada; //pega as coordenadas do ponto
			v2 = vertices[f2].coordenada;
			v3 = vertices[f1].coordenada;


			objetos.push_back(new Triangle(v1, v2, v3, materialTriangle)); //adiciona o objeto ao vetor
		}

	}
	else cout << "Unable to open file";


	myfile.close();
	return objetos;
}

vec rayCast(Ray r, std::vector <Object*> objetos, std::vector<Luz> luzes, int salto) {
	vec cor;
	vec corResto;
	cor << 0.0 << 0.0 << 0.0;
	corResto << 127.0 << 127.0 << 127.0;

	if (salto < 3) {
		Ambient a;
		a.ka << 127.0 << 127.0 << 127.0;

		//luz ambiente
		double intensidade = 0.15;
		vec posicao_luz;
		posicao_luz << 0.0 << 0.0 << -2.0;
		vec radiancia_luz;
		radiancia_luz << intensidade << intensidade << intensidade;
		Luz luzA;
		luzA.pos = posicao_luz;
		luzA.radiancia = radiancia_luz;

		bool intersected = false;
		Intersection da;
		Intersection shadow;

		/*
		l = pos_luz - p;
		Intersection shadow;
		Ray r2(da.p,l);
		for(int o = 0;o<objetos.size();++o){
		if (!(objetos[o]->intersect(r2,&shadow) && shadow.t>0.0 && shadow.t <1.0;
		cor+= blinn phong
		}
		*/

		for (int i = 0; i < objetos.size(); ++i) {
			Intersection da_temp;
			if (objetos[i]->intersect(r, da_temp)) {
				if (intersected && da_temp.t < da.t) {
					da = da_temp;
				}
				else if (!intersected) {
					intersected = true;
					da = da_temp;
				}
			}
		}

		if (intersected) {
			vec k;
			k << 0.0 << 0.0 << 0.0;

			//Intersection reflect;
			//rLinhaDir /= norm(rLinhaDir);
			//reflect.normal /= norm(reflect.normal);

			for (int li = 0; li < luzes.size(); ++li) {
				Luz luz_atual = luzes[li];
				vec l = luz_atual.pos - da.p;
				l /= norm(l);
				da.normal /= norm(da.normal);
				// if (dot(l, da.normal) < 0) da.normal *= -1;  // <<serve para o fix das faces dos triangulos do cubo
				r.direction /= norm(r.direction);
				vec h = (-r.direction + l);
				h /= norm(h);

				Ray r2(da.p, l);

				for (int o = 0; o < objetos.size(); ++o) {
					if (!(objetos[o]->intersect(r2, shadow) && shadow.t > 0.0 && shadow.t < 1.0)) {
						k += ((da.material.kd % luz_atual.radiancia)*std::max(0.0, dot(da.normal, l)) +
							(da.material.ke % luz_atual.radiancia)*pow(std::max(0.0, dot(da.normal, h)), da.material.shininess));
					}
				}
			}

			vec rLinhaDir = (r.direction) - (2.0 * (dot(r.direction, da.normal))*da.normal);
			Ray rLinha(da.p + da.normal*0.001, -rLinhaDir);
			if (salto < 3) {
				++salto;
				cor = (1.0 - da.material.reflection)*k + (da.material.reflection*(rayCast(rLinha, objetos, luzes, salto)));
			}
			else {
				cor = k;
			}
		}
	}
	return cor;

}


int main() {

	ofstream output;
	output.open("imagem.pgm");
	output << "P3" << endl;
	output << "800 600" << endl;
	output << "255" << endl;

	vec origin;
	origin << 0.0 << 0.0 << -3.0;

	vec centro;
	centro << 0.0 << 0.0 << 10.0;
	vec centro2;
	centro2 << 2.0 << 0.0 << 10.0;
	vec centro3;
	centro3 << 1.0 << 1.75 << 10.0;


	//luz ambiente
	double intensidade1 = 0.15;
	vec posicao_luz1;
	posicao_luz1 << 0.0 << 0.0 << -2.0;
	vec radiancia_luz1;
	radiancia_luz1 << intensidade1 << intensidade1 << intensidade1;
	Luz luz1;
	luz1.pos = posicao_luz1;
	luz1.radiancia = radiancia_luz1;

	double intensidade2 = 0.05;
	vec posicao_luz2;
	posicao_luz2 << 0.0 << 4.0 << 3.0;
	vec radiancia_luz2;
	radiancia_luz2 << intensidade2 << intensidade2 << intensidade2;
	Luz luz2;
	luz2.pos = posicao_luz2;
	luz2.radiancia = radiancia_luz2;

	double intensidade3 = 0.25;
	vec posicao_luz3;
	posicao_luz3 << 6.0 << 1.0 << 3.0;
	vec radiancia_luz3;
	radiancia_luz3 << intensidade3 << intensidade3 << intensidade3;
	Luz luz3;
	luz3.pos = posicao_luz3;
	luz3.radiancia = radiancia_luz3;

	double intensidade4 = 0.4;
	vec posicao_luz4;
	posicao_luz4 << -6.0 << 1.0 << 0.0;
	vec radiancia_luz4;
	radiancia_luz4 << intensidade4 << intensidade4 << intensidade4;
	Luz luz4;
	luz4.pos = posicao_luz4;
	luz4.radiancia = radiancia_luz4;

	std::vector <Object*> objetos;
	std::vector<Luz> luzes;

	luzes.push_back(luz1);
	luzes.push_back(luz2);
	luzes.push_back(luz3);
	luzes.push_back(luz4);

	Material material;
	material.kd << 255.0 << 255.0 << 0.0;
	material.ke << 255.0 << 255.0 << 255.0;
	material.shininess = 50.0;
	material.reflection = 0.8;

	Material material2;
	material2.kd << 255.0 << 0.0 << 0.0;
	material2.ke << 255.0 << 255.0 << 255.0;
	material2.shininess = 50.0;
	material2.reflection = 0.5;

	Material materialTriangle;
	materialTriangle.kd << 0.0 << 255.0 << 0.0;
	materialTriangle.ke << 255.0 << 255.0 << 255.0;
	materialTriangle.reflection = 0.9;
	materialTriangle.shininess = 35.0;

	vec a1, a2, a3, a4;
	a1 << -8.0 << -1.0 << 0.0;
	a2 << 8.0 << -1.0 << 0.0;
	a3 << -8.0 << -1.0 << 20.0;
	a4 << 8.0 << -1.0 << 20.0;

	Triangle t1(a1, a2, a3, materialTriangle);
	Triangle t2(a2, a3, a4, materialTriangle);


	objetos.push_back(new Sphere(centro, 1.0, material));
	objetos.push_back(new Sphere(centro2, 1.0, material2));
	objetos.push_back(new Sphere(centro3, 1.0, materialTriangle));
	//objetos.push_back(&t1);
	//objetos.push_back(&t2);


	/*
	vec posicao_luz;
	posicao_luz << 0.0 << 0.0 << -2.0;
	vec radiancia_luz;
	radiancia_luz << intensidade << intensidade << intensidade;
	Luz luzA;
	luzA.pos = posicao_luz;
	luzA.radiancia = radiancia_luz;


	Ambient a;
	a.ka << 127.0 << 127.0 << 127.0;
	*/

	mat k;
	k << 1000.0 << 0.0 << 400.0 << endr
		<< 0.0 << -1000.0 << 300.0 << endr
		<< 0.0 << 0.0 << 1.0;
	mat invk = k.i();

	vec cor;
	cor << 0.0 << 0.0 << 0.0;

	//monta a janela da imagem
	for (int linha = 0; linha < 600; ++linha) {
		for (int coluna = 0; coluna < 800; ++coluna) {
			vec dir;
			dir << coluna << linha << 1.0;
			dir = invk*dir;
			dir = (dir / dir(2))*5.0;
			Ray r(origin, dir);

			//calcula a cor de cada ponto
			cor = rayCast(r, objetos, luzes, 0);

			//imprime a img
			output << std::min(255, int(cor(0))) << " "
				<< std::min(255, int(cor(1))) << " "
				<< std::min(255, int(cor(2))) << " ";
		}
	}


	output.close();
	return 0;
}
