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

	Material& operator = (const Material& mat) {
		ke = mat.ke;
		kd = mat.kd;
		shininess = mat.shininess;
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
	mat rot;
	double sent = 0.0174; //seno de 1 grau, se for grande nao aparece na imagem
	double cost = 0.9998; //cos de 1 grau

	/*rot << cost << 0.0 << -sent << endr //matriz de rotacao
		<< 0.0 << 1.0 << 0.0 << endr
		<< sent << 0.0 << cost;*/


	std::fstream myfile(nomeArq, std::fstream::in);
	if (myfile.is_open()) {
		std::string line;
		while (getline(myfile, line)) {
			std::stringstream stream(line, std::stringstream::in);
			std::string type;
			stream >> type;

			if (type == "#" || type == "") {
				continue;
			}else 
				if (type == "v") {
				float x, y, z;
				stream >> x >> y >> z; //pega as coordenadas x,y,z do arquivo
				v.coordenada << x << y << (z+10); //cria o vertice com as coordenadas
				//v.coordenada = rot * v.coordenada; //multiplica as coordenadas pela matriz de rotacao
				vertices.push_back(v); //adiciona o vertice ao vector

			}else
				if (type == "vn") {
					float t, r, s;
					stream >> t >> r >> s;
					continue;

				}else
					if (type == "f") {
						
						/*int numDeBarras = std::count(line.begin(), line.end(), '/');
						char c = '/';
						cout  << line.find(c) << numDeBarras;
						cin.get();

						if (numDeBarras == 0) {
							stream >> f >> g >> h;
						}
						else if (numDeBarras == 3) {
							int barraN1 = findNext(line, 0, '/' );
							int barraN2 = findNext(line, barraN1, '/');
							f = stoi(line.substr(0, barraN1)); //pega a o primeiro vertice da face
							g = stoi(line.substr((barraN1+1), barraN2));
							h = stoi(line.substr((barraN2), line.length()));
						//Ta tudo errado

						}
						else if (numDeBarras == 6) {

						}*/


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

			cout << "f1 - " << f1 << endl;
			cout << "f2 - " << f2 << endl;
			cout << "f3 - " << f3 << endl;

			cout << "v1 - " << v1 << endl;
			cout << "v2 - " << v2 << endl;
			cout << "v3 - " << v3 << endl;


			objetos.push_back(new Triangle(v1, v2, v3, materialTriangle)); //adiciona o objeto ao vetor
		}

	}
	else cout << "Unable to open file";


	myfile.close();
	//cin.get();
	return objetos;
}


int main() {

	double intensidade = 0.9;
	vec posicao_luz;
	posicao_luz << -5.0 << 4.0 << 8.0;
	vec radiancia_luz;
	radiancia_luz << intensidade << intensidade << intensidade;
	Luz luz;
	luz.pos = posicao_luz;
	luz.radiancia = radiancia_luz;

	mat k;
	k << 1000.0 << 0.0 << 400.0 << endr
		<< 0.0 << -1000.0 << 300.0 << endr
		<< 0.0 << 0.0 << 1.0;
	mat invk = k.i();

	Material materialTriangle;
	materialTriangle.kd << 255.0 << 255.0 << 255.0;
	materialTriangle.ke << 0.0 << 0.0 << 0.0;
	materialTriangle.shininess = 20.0;

	ofstream output;
	output.open("imagem.pgm");
	output << "P3" << endl;
	output << "800 600" << endl;
	output << "255" << endl;

	Material materialSphere;
	materialSphere.kd << 0.0 << 0.0 << 255.0;
	materialSphere.ke << 255.0 << 255.0 << 255.0;
	materialSphere.shininess = 60.0;

	Material materialSphereB;
	materialSphereB.kd << 255.0 << 0.0 << 0.0;
	materialSphereB.ke << 255.0 << 255.0 << 255.0;
	materialSphereB.shininess = 20.0;

	
	vec origin;
	origin << 0.0 << 0.0 << 0.0;
	vec centro;
	centro << 0.0 << 0.0 << 10.0;
	vec centro2;
	centro2 << 2.0 << 0.0 << 10.0;
	
	std::vector<Object*> objetos;
	//objetos = loadObject("cubonovo.txt", materialTriangle); //carregar objeto
	objetos.push_back(new Sphere(centro, 1.0, materialSphere));
	objetos.push_back(new Sphere(centro2, 1.0, materialSphereB));
	vec a1, a2, a3, a4;
	a1 << -8.0 << -1.0 << 0.0;
	a2 << 8.0 << -1.0 << 0.0;
	a3 << -8.0 << -1.0 << 20.0;
	a4 << 8.0 << -1.0 << 20.0;

	Triangle t1(a1, a2, a3, materialSphere);
	Triangle t2(a2, a3, a4, materialSphere);
	std::vector<Luz> luzes;
	luzes.push_back(luz);
	objetos.push_back(&t1);
	objetos.push_back(&t2);

	for (int linha = 0; linha < 600; ++linha) {
		for (int coluna = 0; coluna < 800; ++coluna) {
			bool intersected = false;
			vec dir;
			dir << coluna << linha << 1.0;
			dir = invk*dir;
			dir = (dir / dir(2))*5.0;
			Ray r(origin, dir);
			Intersection da;
			for (int i = 0; i<objetos.size(); ++i) {
				Intersection da_temp;
				if (objetos[i]->intersect(r, da_temp)) {
					if (intersected && da_temp.t<da.t) {
						da = da_temp;
					}
					else if (!intersected) {
						intersected = true;
						da = da_temp;
					}
				}
			}

			if (intersected) {
				vec cor;
				cor << 0.0 << 0.0 << 0.0;
				for (int li = 0; li < luzes.size(); ++li) {
					Luz luz_atual = luzes[li];
					vec l = luz_atual.pos - da.p;
					
					da.normal /= norm(da.normal);
					if (da.flipNormal && dot(da.normal, l) < 0.0) da.normal *= -1.0;
					
					Intersection shadow;
					Ray r2(da.p, l);
					l /= norm(l);
					vec h = l + da.normal;
					h /= norm(h);
					bool paint = true;
					for (int o = 0; o < objetos.size(); ++o) {
						bool hit = objetos[o]->intersect(r2, shadow);
						if (hit && (shadow.t > 0.0) && (shadow.t < 1.0)) {
							paint = false;
							break;
						}
					}
					if (paint) {
						cor += (da.material.kd % luz_atual.radiancia)*std::max(0.0, dot(da.normal, l)) +
							(da.material.ke % luz_atual.radiancia)*pow(std::max(0.0, dot(h, da.normal)), da.material.shininess);
					}
					
				}

				output << std::min(255, int(cor(0))) << " "
					<< std::min(255, int(cor(1))) << " "
					<< std::min(255, int(cor(2))) << " ";
			}
			else {
				output << "255 255 255 ";
			}
		}
	}
	output.close();
	return 0;
}
