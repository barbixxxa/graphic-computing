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
	Intersection(double t, const vec& p, const vec& normal, Material material) : t(t), p(p), normal(normal), material(material) {}
	Intersection() : t(-1.0) {}

	Intersection& operator = (const Intersection& i) {
		t = i.t;
		p = i.p;
		normal = i.normal;
		material = i.material;
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

		//t
		/*M4.insert_cols(0, A_ - B_);
		M4.insert_cols(1, A_ - C_);
		M4.insert_cols(2, A_ - r.origin);*/

		double detM = det(M);
		double beta = det(M2) / (detM);
		double gamma = det(M3) / (detM);
		//da.t = det(M4) / (detM);

		const vec& ba = B_ - A_;
		const vec& ca = C_ - A_;
		const vec& oa = r.origin - A_;
		da.normal = cross(ba, ca);
		da.t = (-(dot(da.normal, oa))) / (dot(r.direction, da.normal));


		if ((gamma<0) || (gamma>1)) return false;
		if ((beta<0) || (beta>(1.0 - gamma))) return false;

		//da.normal = cross(B_ - A_, C_ - A_);
		da.p = r.origin + r.direction *da.t;
		da.material = material_;
		return true;
	}
};


int main() {
	double intensidade = 0.9;
	vec posicao_luz;
	posicao_luz << -5.0 << 8.0 << 0.0;
	vec radiancia_luz;
	radiancia_luz << intensidade << intensidade << intensidade;
	Luz luz;
	luz.pos = posicao_luz;
	luz.radiancia = radiancia_luz;

	Luz luz2;
	luz2.pos << 5.0 << -8.0 << 0.0;
	luz2.radiancia << 0.1 << 0.9 << 0.1;

	Luz luz3;
	luz3.pos << 5.0 << -8.0 << 2.0;
	luz3.radiancia << 0.9 << 0.1 << 0.1;

	Luz luz4;
	luz4.pos << 5.0 << -8.0 << 2.0;
	luz4.radiancia << 0.1 << 0.1 << 0.9;

	mat k;
	k << 1000.0 << 0.0 << 400.0 << endr
		<< 0.0 << -1000.0 << 300.0 << endr
		<< 0.0 << 0.0 << 1.0;
	mat invk = k.i();

	Material materialSphere;
	materialSphere.kd << 0.0 << 0.0 << 255.0;
	materialSphere.ke << 255.0 << 255.0 << 255.0;
	materialSphere.shininess = 60.0;

	Material materialTriangle;
	materialTriangle.kd << 255.0 << 0.0 << 0.0;
	materialTriangle.ke << 0.0 << 0.0 << 0.0;
	materialTriangle.shininess = 60.0;

	Material materialTriangleD;
	materialTriangleD.kd << 0.0 << 255.0 << 0.0;
	materialTriangleD.ke << 0.0 << 0.0 << 0.0;
	materialTriangleD.shininess = 60.0;

	Material materialTriangleT;
	materialTriangleT.kd << 0.0 << 0.0 << 0.0;
	materialTriangleT.ke << 0.0 << 0.0 << 0.0;
	materialTriangleT.shininess = 60.0;

	ofstream output;
	output.open("imagem.pgm");
	output << "P3" << endl;
	output << "800 600" << endl;
	output << "255" << endl;

	vec origin;
	origin << 0.0 << 0.0 << 0.0;
	vec centro;
	centro << 0.0 << 0.0 << 10.0;
	vec A, B, C, D, E, F, G, H, I;
	A << 4.0 << 0.0 << 10.0;
	B << 0.0 << 3.0 << 5.0;
	C << 0.0 << 0.0 << 10.0;

	D << -4.0 << 0.0 << 10.0;
	E << 0.0 << -3.0 << 10.0;
	F << 0.0 << 0.0 << 10.0;

	G << -4.0 << -2.0 << 10.0;
	H << 0.0 << 3.0 << 5.0;
	I << 0.0 << 0.0 << 10.0;

	std::vector <Object*> objetos;
	std::vector<Luz> luzes;
	objetos.push_back(new Sphere(centro, 2.0, materialSphere));
	objetos.push_back(new Triangle(A, B, C, materialTriangle));
	objetos.push_back(new Triangle(D, E, F, materialTriangleD));
	objetos.push_back(new Triangle(G, H, I, materialTriangleT));
	luzes.push_back(luz);
	luzes.push_back(luz2);
	luzes.push_back(luz3);
	luzes.push_back(luz4);

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
					l /= norm(l);
					da.normal /= norm(da.normal);
					if (dot(da.normal, l) < 0.0) da.normal *= -1.0;
					vec h = (luz_atual.pos - da.p) + da.normal;
					h /= norm(h);
					cor += (da.material.kd % luz_atual.radiancia)*std::max(0.0, dot(da.normal, l)) +
						(da.material.ke % luz_atual.radiancia)*pow(std::max(0.0, dot(h, da.normal)), da.material.shininess);
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



// Objetivo = leitor de .obj (arquivo vindo do blender,maya...)
