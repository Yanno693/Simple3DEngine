#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
//#include <iostream>

// A 3D Vector with (x,y,z) coordinates
class Vec3
{
	private:
		double x, y, z; // coordinates

	public:
		Vec3()
		{
			x = 0;
			y = 0;
			z = 0;
		}

		Vec3(double _x, double _y, double _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}

		double getX() const { return x; }

		double getY() const { return y; }

		double getZ() const { return z; }

		void setX(const double _x) { x = _x; }

		void setY(const double _y) { y = _y; }

		void setZ(const double _z) { z = _z; }

		Vec3 operator+(const Vec3& v)
		{
			return Vec3(x + v.getX(), y + v.getY(), z + v.getZ());
		}

		Vec3& operator+=(const Vec3& v)
		{
			*this = *this + v;
			return *this;
		}

		Vec3 operator-(const Vec3& v)
		{
			return Vec3(x - v.getX(), y - v.getY(), z - v.getZ());
		}

		Vec3& operator-=(const Vec3& v)
		{
			*this = *this - v;
			return *this;
		}

		Vec3 operator*(const double& n)
		{
			return Vec3(x*n, y*n, z*n);
		}

		Vec3& operator*=(const double& n)
		{
			*this = *this * n;
			return *this;
		}

		Vec3 operator/(const double& n)
		{
			if (n == 0.0)
			{
				double nSafe = n + 0.001;
				return Vec3(x / nSafe, y / nSafe, z / nSafe);
			}
			else
			{
				return Vec3(x / n, y / n, z / n);
			}
		}

		Vec3& operator/=(const double& n)
		{
			*this = *this / n;
			return *this;
		}

		friend std::ostream& operator<< (std::ostream& stream, Vec3& v)
		{
			return stream << "(" << v.x << "," << v.y << "," << v.z << ")";
		}
};

// A triangle defined by three points
class Triangle
{
	private:
		Vec3 p[3]; // vertices

	public:
		Triangle()
		{
			p[0] = Vec3(0, 0, 0);
			p[1] = Vec3(0, 1, 0);
			p[2] = Vec3(1, 1, 0);
		}

		Triangle(Vec3 p0, Vec3 p1, Vec3 p2)
		{
			p[0] = p0;
			p[1] = p1;
			p[2] = p2;
		}

		Vec3 get(int i) const
		{
			return p[i];
		}

		void set(int i, Vec3 v)
		{
			p[i] = v;
		}

		Vec3& operator[](const int i) { return p[i]; }

		friend std::ostream& operator<< (std::ostream& stream, Triangle& t)
		{
			return stream << '<' << t.p[0] << ',' << t.p[1] << ',' << t.p[2] << '>';
		}
};

// A Mesh made of triangles
class Mesh
{
	private:

		std::vector<Triangle> triangles;
		Vec3 position;
		Vec3 scale;

	public:

		Mesh()
		{
			position = Vec3(0, 0, 0);
			scale = Vec3(1, 1, 1);
		}

		void add(const Triangle& tr) { triangles.push_back(tr); }

		void setPosition(const Vec3& _pos) { position = _pos; }

		Vec3 getPosition() const { return position; }

		void translate(const Vec3& _pos) { position += _pos; }

		void setScale(const Vec3& _scale) { scale = _scale; }

		Vec3 getScale() const { return scale; }

		Triangle getTriangle(int i) const { return triangles[i]; }

		Triangle& operator[](const unsigned int i) { return triangles[i]; }

		unsigned int size() const { return triangles.size(); }

		friend std::ostream& operator<< (std::ostream& stream, Mesh& m)
		{
			return stream << "Mesh : [" << m.position << "," << m.scale << "]";
		}
};

// A 4x4 Matrix (for projection of a point on the screen)
class Matrix4x4
{
	private:
		double* m[4];

	public:
		Matrix4x4()
		{
			for (int i = 0; i < 4; i++)
			{
				m[i] = new double[4];
				for (int j = 0; j < 4; j++)
					m[i][j] = 0.0;
			}
		}

		~Matrix4x4()
		{
			for (int i = 0; i < 4; i++)
				delete[] m[i];
		}

		double* operator[](const unsigned int i)
		{
			return m[i];
		}

		Vec3 mul(const Vec3& v) const
		{
			Vec3 res;
			res.setX(v.getX() * m[0][0] + v.getY() * m[1][0] + v.getZ() * m[2][0] + m[3][0]);
			res.setY(v.getX() * m[0][1] + v.getY() * m[1][1] + v.getZ() * m[2][1] + m[3][1]);
			res.setZ(v.getX() * m[0][2] + v.getY() * m[1][2] + v.getZ() * m[2][2] + m[3][2]);

			double w = v.getX() * m[0][3] + v.getY() * m[1][3] + v.getZ() * m[2][3] + m[3][3];

			if (w != 0.0)
			{
				res /= w;
			}

			return res;
		}

		friend std::ostream& operator<< (std::ostream& stream, Matrix4x4& mx)
		{
			for (int i = 0; i < 4; i++)
				stream << '[' << mx.m[i][0] << ',' << mx.m[i][1] << ',' << mx.m[i][2] << ',' << mx.m[i][3] << ']' << (i != 3 ? "\n" : " ");

			return stream;
		}
};

class GameEngine : public olc::PixelGameEngine
{
	private:
		Matrix4x4 projection;
		float delta;
		Vec3 camera;

		std::vector<Mesh> mesh;

		void show(const Mesh& m)
		{
			for (unsigned int i = 0; i < m.size(); i++)
			{
				Triangle tScale;
				for (unsigned int t = 0; t < 3; t++)
				{
					tScale[t].setX(m.getTriangle(i)[t].getX() * m.getScale().getX());
					tScale[t].setY(m.getTriangle(i)[t].getY() * m.getScale().getY());
					tScale[t].setZ(m.getTriangle(i)[t].getZ() * m.getScale().getZ());
				}
				
				Triangle tPosition = tScale;
				for (unsigned int t = 0; t < 3; t++)
				{
					tPosition[t] += m.getPosition() - camera;
				}

				Triangle tProjection = Triangle(
					projection.mul(tPosition[0]),
					projection.mul(tPosition[1]),
					projection.mul(tPosition[2]) );

				for (unsigned int t = 0; t < 3; t++)
				{
					tProjection[t] = tProjection[t] + Vec3(1, 1, 0);

					tProjection[t].setX(tProjection[t].getX() * 0.5 * ScreenWidth());
					tProjection[t].setY(tProjection[t].getY() * 0.5 * ScreenHeight()); //+= Vec3(0.5 * ScreenWidth(), 0.5 * ScreenHeight(), 0);
				}

				DrawTriangle(
					tProjection[0].getX(),
					tProjection[0].getY(),
					tProjection[1].getX(),
					tProjection[1].getY(),
					tProjection[2].getX(),
					tProjection[2].getY(),
					olc::BLACK
				);
			}
		}

	public : 
		
		GameEngine()
		{
			sAppName = "Easy Fast Game Engine";
			delta = 0.0;
			camera = Vec3(0, 0, 0);
		}

		void addMesh(const Mesh& m)
		{
			mesh.push_back(m);
		}

		bool OnUserCreate() override
		{
			double fNear = 0.1;
			double fFar = 1000.0;
			double fov = 90.0;

			double ar = ScreenWidth() / ScreenHeight();
			double frad = 1.0 / tan(fov * 0.5 / 180.0 * 3.14159);

			projection[0][0] = ar * frad;
			projection[1][1] = frad;
			projection[2][2] = fFar / (fFar - fNear);
			projection[3][2] = (-fFar * fNear) / (fFar - fNear);
			projection[2][3] = 1.0;
			
			Mesh m = Mesh();
			
			m.add(Triangle(Vec3(0, 0, 0), Vec3(0, 1, 0), Vec3(1, 1, 0)));
			m.add(Triangle(Vec3(0, 0, 0), Vec3(1, 1, 0), Vec3(1, 0, 0)));

			m.add(Triangle(Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(1, 1, 1)));
			m.add(Triangle(Vec3(1, 0, 0), Vec3(1, 1, 1), Vec3(1, 0, 1)));

			m.add(Triangle(Vec3(1, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1)));
			m.add(Triangle(Vec3(1, 0, 1), Vec3(0, 1, 1), Vec3(0, 0, 1)));

			m.add(Triangle(Vec3(0, 0, 1), Vec3(0, 1, 1), Vec3(0, 1, 0)));
			m.add(Triangle(Vec3(0, 0, 1), Vec3(0, 1, 0), Vec3(0, 0, 0)));

			m.add(Triangle(Vec3(0, 1, 0), Vec3(0, 1, 1), Vec3(1, 1, 1)));
			m.add(Triangle(Vec3(0, 1, 0), Vec3(1, 1, 1), Vec3(1, 1, 0)));

			m.add(Triangle(Vec3(1, 0, 1), Vec3(0, 0, 1), Vec3(0, 0, 0)));
			m.add(Triangle(Vec3(1, 0, 1), Vec3(0, 0, 0), Vec3(1, 0, 0)));
			m.setPosition(Vec3(0,0,3));

			Mesh m2 = m;

			m2.translate(Vec3(1, 1, 1));

			Mesh m3 = m;

			m3.translate(Vec3(2, 2, 2));
			m3.setScale(Vec3(2, 2, 2));

			addMesh(m);
			addMesh(m2);
			addMesh(m3);

			return true;
		}

		bool OnUserUpdate(float elapsedTime) override
		{
			FillRect(0, 0, ScreenWidth(), ScreenHeight(), olc::WHITE);

			delta += 1.0f * elapsedTime;

			if (GetKey(olc::Key::UP).bPressed)
			{
				camera += Vec3(0, 0, 1);
			}

			if (GetKey(olc::Key::DOWN).bPressed)
			{
				camera += Vec3(0, 0, -1);
			}

			if (GetKey(olc::Key::LEFT).bPressed)
			{
				camera += Vec3(-1, 0, 0);
			}

			if (GetKey(olc::Key::RIGHT).bPressed)
			{
				camera += Vec3(1, 0, 0);
			}

			if (GetKey(olc::Key::SPACE).bPressed)
			{
				camera += Vec3(0, -1, 0);
			}

			if (GetKey(olc::Key::CTRL).bPressed)
			{
				camera += Vec3(0, 1, 0);
			}

			for (unsigned int i = 0; i < mesh.size(); i++)
			{
				show(mesh[i]);
			}

			/*std::vector<std::thread> vt;

			for (unsigned int i = 0; i < mesh.size(); i++) // Not faster ? :/
			{
				vt.push_back(std::thread(&GameEngine::show, this, mesh[i]));//
			}

			for (unsigned int i = 0; i < mesh.size(); i++)
			{
				vt[i].join();
			}*/

			return true;
		}
};

int main()
{
	using namespace std;

	GameEngine ge;
	if (ge.Construct(200, 200, 4, 4))
		ge.Start();
	else
		cout << "Error" << endl;

	return 0;
}

/*int main()
{
	using namespace std;
	
	Vec3 v;

	cout << v << endl;

	v += Vec3(1, 2, 3);

	cout << v << endl;

	v *= 2;

	cout << v << endl;

	v -= Vec3(1,1,1);

	cout << v << endl;

	Triangle t;

	cout << t[0] << endl;

	t[0] += Vec3(1,1,1);

	cout << t[0] << endl;

	Triangle t2 = t;

	cout << t2[0] << endl;

	t[0] += Vec3(1, 1, 1);

	cout << t2[0] << endl;
	cout << t[0] << endl;

	Vec3 vv = t[0];

	cout << t[0] << " " << vv << endl;
	t[0] += Vec3(1, 1, 1);

	cout << t[0] << " " << vv << endl;

	cout << t << endl;

	cout << t2 << endl;

	return 0;
}*/