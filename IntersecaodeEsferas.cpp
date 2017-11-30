#include <fstream>
#include <armadillo>
#include <thread>

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

int main() {
	vec centro_esfera;
	centro_esfera <<  0.0 << 0.0 << 10.0;
	vec centro_esfera2;
	centro_esfera2 << 2.0 << 0.0 << 10.0;
	vec centro_esfera3;
	centro_esfera3 << -2.0 << 0.0 << 10.0;
	vec origin;
	origin << 0.0 << 0.0 << 0.0;
	double intensidade = 0.9;
	vec posicao_luz;
	posicao_luz << 0.0 << 4.0 << 4.0 ;
	vec radiancia_luz;
	radiancia_luz << intensidade << intensidade << intensidade ;
	mat k;
	k << 1000.0 << 0.0 << 400.0 << endr
		<< 0.0 << -1000.0 << 300.0 << endr
		<< 0.0 << 0.0 << 1.0;
	mat invk = k.i();

	ofstream output;
	output.open("imagem.pgm");
	output << "P3" << endl;
	output << "800 600" << endl;
	output << "255" << endl;

	for (int linha = 0; linha < 600; linha++) {
		for (int coluna = 0; coluna < 800; coluna++) {
			vec r;
			r << double(coluna) << double(linha) << 1.0 ;
			
			vec normal, hit;

			r = invk * r;
			r *= 5.0 / r(2);
			if (intersect(origin, r, centro_esfera, 2.0, normal, hit)) {
				normal /= norm(normal);
				vec l = posicao_luz - hit;
				l /= norm(l);

				vec cor;
				cor << 255.0 << 0.0 << 0.0;
				cor = (cor % radiancia_luz)*std::max(0.0, dot(normal, l));
				output << std::min(255, int(cor(0))) << " "
					<< std::min(255, int(cor(1))) << " "
					<< std::min(255, int(cor(2))) << " ";
			}else if (intersect(origin, r, centro_esfera2, 1.0, normal, hit)) {
				normal /= norm(normal);
				vec l = posicao_luz - hit;
				l /= norm(l);

				vec cor;
				cor << 0.0 << 0.0 << 255.0;
				cor = (cor % radiancia_luz)*std::max(0.0, dot(normal, l));
				output << std::min(255, int(cor(0))) << " "
					<< std::min(255, int(cor(1))) << " "
					<< std::min(255, int(cor(2))) << " ";
			}else if (intersect(origin, r, centro_esfera3, 1.0, normal, hit)) {
				normal /= norm(normal);
				vec l = posicao_luz - hit;
				l /= norm(l);

				vec cor;
				cor << 0.0 << 255.0 << 0.0;
				cor = (cor % radiancia_luz)*std::max(0.0, dot(normal, l));
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

class Object {
private:
	vec kd_;
	vec ke_;
	double shininess_;

public:
	

	vec kd() const {
		return kd_;
	}

	vec ke() const {
		return ke_;
	}

	double shininess() const {
		return shininess_;
	}

	struct Ray {
		vec origin;
		vec direction;
	};

	struct Intersection {
		vec normal;
		double t;
		vec p;
	};
	virtual bool intersect(const Ray& r, Intersection& da) const = 0;
};

struct Material {
	vec kd;
	vec ke;
	double shininess;
};

struct Luz {
	vec pos;
	vec radiancia;
};

struct Intersection {
	double t;
	vec p;
	vec normal;
	Material material;
};

struct Ray {
	vec origin;
	vec direction;
};

};

class Sphere : public Object {
private:
	vec center_;
	double radius_;
public:
	Sphere() {

	}

	virtual bool intersect(const Ray& r, Intersection& da) const {

		vec oc = r.origin - center_;
		double a = dot(r.direction, r.direction);
		double b = 2.0*dot(r.direction, oc);
		double c = dot(oc, oc) - radius_*radius_;

		double delta = b*b - 4.0*a*c;

		if (delta < 0.0) { return false; }

		double t1 = (-b + sqrt(delta) / 2.0*a);
		double t2 = (-b - sqrt(delta) / 2.0*a);

		da.t = (t1 < t2) ? t1 : t2;
		da.p = r.origin + r.direction *da.t;
		da.normal = da.p - center_;

		return true;
	}
};

class Triangle : public Object {
private:
	vec A, B, C;
public:
	Triangle() {

	}

	virtual bool intersect(const Ray& r, Intersection& da) const {

		mat M, M2, M3, M4;

		//A
		M.insert_cols(0, A - B);
		M.insert_cols(1, A - C);
		M.insert_cols(2, r.direction);

		//beta
		M2.insert_cols(0, A - r.origin);
		M2.insert_cols(1, A - C);
		M2.insert_cols(2, r.direction);

		//gamma
		M3.insert_cols(0, A - B);
		M3.insert_cols(1, A - r.origin);
		M3.insert_cols(2, r.direction);

		//t
		M4.insert_cols(0, A - B);
		M4.insert_cols(1, A - C);
		M4.insert_cols(2, A - r.origin);

		double detM = det(M);
		double beta = det(M2) / detM;
		double gamma = det(M3) / detM;
		da.t = det(M4) / detM;

		if ((gamma<0) || (gamma>1)) return false;
		if ((beta<0) || (beta>(1.0 - gamma))) return false;

		da.normal = cross(B - A, C - A);
		da.p = r.origin + r.direction *da.t;

		return true;
	}


};






// Objetivo = leitor de .obj (arquivo vindo do blender,maya...)
