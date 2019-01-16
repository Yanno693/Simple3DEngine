#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <algorithm>
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

		double magnitude() const
		{
			return sqrt(x*x + y * y + z * z);
		}

		void normalize()
		{
			double w = magnitude();
			
			if(w != 0.0)
				(*this) /= w;
		}

		void rotateX(const double angle) // angle in rad
		{
			double _y, _z;
			_y = y;
			_z = z;
			
			y = _y * cos(angle) + _z * sin(angle);
			z = _y * (-sin(angle)) + _z * cos(angle);
		}

		void rotateY(const double angle) // angle in rad
		{
			double _x, _z;
			_x = x;
			_z = z;

			x = _x * cos(angle) + _z * (-sin(angle));
			z = _x * sin(angle) + _z * cos(angle);
		}

		void rotateZ(const double angle) // angle in rad
		{
			double _x, _y;
			_x = x;
			_y = y;

			x = _x * cos(angle) + _y * sin(angle);
			y = _x * (-sin(angle)) + _y * cos(angle);
		}

		double dot(const Vec3& v)
		{
			return (x * v.getX() + y * v.getY() + z * v.getZ());
		}

		friend std::ostream& operator<< (std::ostream& stream, Vec3 v)
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

		Vec3 normal()
		{
			Vec3 res, v1, v2;

			v1 = p[1] - p[0];
			v2 = p[2] - p[0];
			
			res.setX(v1.getY() * v2.getZ() - v1.getZ() * v2.getY());
			res.setY(v1.getZ() * v2.getX() - v1.getX() * v2.getZ());
			res.setZ(v1.getX() * v2.getY() - v1.getY() * v2.getX());

			res.normalize();

			return res;
		}

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

		// open mesh from Wavefront .obj file (with triangulated faces) (https://en.wikipedia.org/wiki/Wavefront_.obj_file) 
		Mesh(std::string src)
		{
			position = Vec3(0, 0, 0);
			scale = Vec3(1, 1, 1);

			std::ifstream file;			
			file.open(src.c_str());
			if (file.is_open())
			{
				std::string s;
				std::vector<Vec3> vP; // Values of each vertices of the mesh
				std::vector<int> vT; // Indexes of the vertices in the previous array for each triangle
				while (!file.eof())
				{
					std::getline(file, s);
					if (!file.eof())
					{						
						if (s[0] == 'v')
						{							
							double x, y, z;
							std::string::size_type c1, c2;

							x = std::stod(s.substr(2),&c1);
							y = std::stod(s.substr(2 + c1), &c2);
							z = std::stod(s.substr(2 + c2 + c2));

							Vec3 p = Vec3(x, y, z);
							vP.push_back(p);
						}

						else if (s[0] == 'f')
						{
							int x, y, z;
							std::string::size_type c1, c2;

							x = std::stoi(s.substr(2), &c1); x -= 1;
							y = std::stoi(s.substr(2 + (c1+1)), &c2); y -= 1;
							z = std::stoi(s.substr(2 + (c1+c2+2))); z -= 1;

							vT.push_back(x);
							vT.push_back(y);
							vT.push_back(z);
						}
					}
				}

				file.close();

				for (unsigned int i = 0; i < vT.size(); i += 3)
				{
					Triangle t = Triangle();
					t[0] = vP[vT[i]];
					t[1] = vP[vT[i + 1]];
					t[2] = vP[vT[i + 2]];

					triangles.push_back(t);
				}
			}
			else
			{
				std::cout << "file could not open" << std::endl;
				exit(1);
			}
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
class Matrix4
{
	private:
		double* m[4];

	public:
		Matrix4()
		{
			for (int i = 0; i < 4; i++)
			{
				m[i] = new double[4];
				for (int j = 0; j < 4; j++)
					m[i][j] = 0.0;
			}
		}
	
		~Matrix4()
		{
			for (int i = 0; i < 4; i++)
				delete[] m[i];
		}

		double* operator[](const unsigned int i) { return m[i]; }

		Vec3 mul(const Vec3& v) const
		{
			Vec3 res;
			res.setX(v.getX() * m[0][0] + v.getY() * m[1][0] + v.getZ() * m[2][0] + m[3][0]);
			res.setY(v.getX() * m[0][1] + v.getY() * m[1][1] + v.getZ() * m[2][1] + m[3][1]);
			res.setZ(v.getX() * m[0][2] + v.getY() * m[1][2] + v.getZ() * m[2][2] + m[3][2]);

			double w = v.getX() * m[0][3] + v.getY() * m[1][3] + v.getZ() * m[2][3] + m[3][3];

			if (w != 0.0)
				res /= w;

			return res;
		}

		friend std::ostream& operator<< (std::ostream& stream, Matrix4& mx)
		{
			for (int i = 0; i < 4; i++)
				stream << '[' << mx.m[i][0] << ',' << mx.m[i][1] << ',' << mx.m[i][2] << ',' << mx.m[i][3] << ']' << (i != 3 ? "\n" : " ");
				return stream;
		}
};

class GameEngine : public olc::PixelGameEngine
{
	private:
		Matrix4 projection;
		float delta;
		Vec3 cameraPosition;
		Vec3 cameraRotation;
		Vec3 light;

		std::vector<Mesh> mesh;

		std::vector<Triangle> distanceBuffer;
		//std::vector<>


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

				/*for (unsigned int t = 0; t < 3; t++)
				{
					tScale[t].rotateZ(delta * 0.5);
				}*/
				
				Triangle tPosition = tScale;
				for (unsigned int t = 0; t < 3; t++)
				{
					tPosition[t] += m.getPosition() - cameraPosition;
				}

				for (unsigned int t = 0; t < 3; t++)
				{
					tPosition[t].rotateX(-cameraRotation.getX());
					tPosition[t].rotateY(-cameraRotation.getY());
					tPosition[t].rotateZ(-cameraRotation.getZ());
				}

				if (tPosition.normal().dot(tPosition[0]) < 0)
				{
					distanceBuffer.push_back(tPosition);
				}
			}
		}
		
		void clipping(Triangle& t, double lightIntensity)
		{
			std::vector<int> outIndexes;
			for (unsigned int v = 0; v < 3; v++)
			{
				if (t[v].getX() < 0 || t[v].getX() > ScreenWidth() || t[v].getY() < 0 || t[v].getY() > ScreenHeight())
					outIndexes.push_back(v);
			}

			std::cout << outIndexes.size() << std::endl;

			// This will be a recursif function i guess, so i need to verify every possibilities
			// If 0, ok
			if (outIndexes.size() == 0)
			{
				FillTriangle(
					t[0].getX(),
					t[0].getY(),
					t[1].getX(),
					t[1].getY(),
					t[2].getX(),
					t[2].getY(),
					olc::Pixel(255 * lightIntensity, 255 * lightIntensity, 255 * lightIntensity)
				);

				DrawTriangle(
					t[0].getX(),
					t[0].getY(),
					t[1].getX(),
					t[1].getY(),
					t[2].getX(),
					t[2].getY(),
					olc::GREEN
				);
			}
			else if (outIndexes.size() == 1)
			{
				Vec3 p1, p2;
				
				if (t[outIndexes[0]].getX() < 0.0) // Left
				{
					Vec3 dirP1 = t[(outIndexes[0] + 1) % 3] - t[outIndexes[0]];
					dirP1 /= dirP1.getY();
					double mulP1 = t[outIndexes[0]].getX() / dirP1.getX();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setX(0.0);

					Vec3 dirP2 = t[(outIndexes[0] + 2) % 3] - t[outIndexes[0]];
					dirP2 /= dirP2.getY();
					double mulP2 = t[outIndexes[0]].getX() / dirP2.getX();
					p2 = t[outIndexes[0]] - (dirP2 * mulP2);
					p2.setX(0.0);

					Triangle t1 = Triangle(t[(outIndexes[0] + 1) % 3], t[(outIndexes[0] + 2) % 3], p1);
					Triangle t2 = Triangle(p2,p1, t[(outIndexes[0] + 2) % 3]);

					clipping(t1, lightIntensity);
					clipping(t2, lightIntensity);
				}
				else if (t[outIndexes[0]].getX() > ScreenWidth()) // Right
				{	
					Vec3 dirP1 = t[(outIndexes[0] + 1) % 3] - t[outIndexes[0]];
					dirP1 /= dirP1.getY();
					double mulP1 = (t[outIndexes[0]].getX() - ScreenWidth()) / dirP1.getX();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setX(ScreenWidth());

					Vec3 dirP2 = t[(outIndexes[0] + 2) % 3] - t[outIndexes[0]];
					dirP2 /= dirP2.getY();
					double mulP2 = (t[outIndexes[0]].getX() - ScreenWidth()) / dirP2.getX();
					p2 = t[outIndexes[0]] - (dirP2 * mulP2);
					p2.setX(ScreenWidth());

					Triangle t1 = Triangle(t[(outIndexes[0] + 1) % 3], t[(outIndexes[0] + 2) % 3], p1);
					Triangle t2 = Triangle(p2, p1, t[(outIndexes[0] + 2) % 3]);

					clipping(t1, lightIntensity);
					clipping(t2, lightIntensity);
				}
				else if (t[outIndexes[0]].getY() < 0.0) // Up
				{
					Vec3 dirP1 = t[(outIndexes[0] + 1) % 3] - t[outIndexes[0]];
					dirP1 /= dirP1.getX();
					double mulP1 = t[outIndexes[0]].getY() / dirP1.getY();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setY(0.0);

					Vec3 dirP2 = t[(outIndexes[0] + 2) % 3] - t[outIndexes[0]];
					dirP2 /= dirP2.getY();
					double mulP2 = t[outIndexes[0]].getY() / dirP2.getY();
					p2 = t[outIndexes[0]] - (dirP2 * mulP2);
					p2.setY(0.0);

					Triangle t1 = Triangle(t[(outIndexes[0] + 1) % 3], t[(outIndexes[0] + 2) % 3], p1);
					Triangle t2 = Triangle(p2, p1, t[(outIndexes[0] + 2) % 3]);

					clipping(t1, lightIntensity);
					clipping(t2, lightIntensity);
				}
				else if (t[outIndexes[0]].getY() > ScreenHeight()) // Down
				{
					Vec3 dirP1 = t[(outIndexes[0] + 1) % 3] - t[outIndexes[0]];
					dirP1 /= dirP1.getX();
					double mulP1 = (t[outIndexes[0]].getY() - ScreenHeight()) / dirP1.getY();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setY(ScreenHeight());

					Vec3 dirP2 = t[(outIndexes[0] + 2) % 3] - t[outIndexes[0]];
					dirP2 /= dirP2.getY();
					double mulP2 = (t[outIndexes[0]].getY() - ScreenHeight()) / dirP2.getY();
					p2 = t[outIndexes[0]] - (dirP2 * mulP2);
					p2.setY(ScreenHeight());

					Triangle t1 = Triangle(t[(outIndexes[0] + 1) % 3], t[(outIndexes[0] + 2) % 3], p1);
					Triangle t2 = Triangle(p2, p1, t[(outIndexes[0] + 2) % 3]);

					clipping(t1, lightIntensity);
					clipping(t2, lightIntensity);
				}
			}
			else if (outIndexes.size() == 2)
			{
				int inIndex = 0;
				if (inIndex == outIndexes[0])
					inIndex++;
				if (inIndex == outIndexes[1])
					inIndex++;

				Vec3 p1, p2;
				
				if(t[outIndexes[0]].getX() < 0.0 && t[outIndexes[1]].getX() < 0.0)
				{
					Vec3 dirP1 = t[inIndex] - t[outIndexes[0]];
					dirP1 /= dirP1.getY();
					double mulP1 = t[outIndexes[0]].getX() / dirP1.getX();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setX(0.0);

					Vec3 dirP2 = t[inIndex] - t[outIndexes[1]];
					dirP2 /= dirP2.getY();
					double mulP2 = t[outIndexes[1]].getX() / dirP2.getX();
					p2 = t[outIndexes[1]] - (dirP2 * mulP2);
					p2.setX(0.0);

					Triangle _t = Triangle(t[inIndex],p2,p1 );
					clipping(_t, lightIntensity);	
				}
				else if (t[outIndexes[0]].getY() < 0.0 && t[outIndexes[1]].getY() < 0.0)
				{
					Vec3 dirP1 = t[inIndex] - t[outIndexes[0]];
					dirP1 /= dirP1.getX();
					double mulP1 = t[outIndexes[0]].getY() / dirP1.getY();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setY(0.0);

					Vec3 dirP2 = t[inIndex] - t[outIndexes[1]];
					dirP2 /= dirP2.getX();
					double mulP2 = t[outIndexes[1]].getY() / dirP2.getY();
					p2 = t[outIndexes[1]] - (dirP2 * mulP2);
					p2.setY(0.0);

					Triangle _t = Triangle(t[inIndex], p2, p1);
					clipping(_t, lightIntensity);
				}
				else if (t[outIndexes[0]].getX() > ScreenWidth() && t[outIndexes[1]].getX() > ScreenWidth())
				{
					Vec3 dirP1 = t[inIndex] - t[outIndexes[0]];
					dirP1 /= dirP1.getY();
					double mulP1 = (t[outIndexes[0]].getX() - ScreenWidth()) / dirP1.getX();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setX(ScreenWidth());

					Vec3 dirP2 = t[inIndex] - t[outIndexes[1]];
					dirP2 /= dirP2.getY();
					double mulP2 = (t[outIndexes[1]].getX() - ScreenWidth()) / dirP2.getX();
					p2 = t[outIndexes[1]] - (dirP2 * mulP2);
					p2.setX(ScreenWidth());

					Triangle _t = Triangle(t[inIndex], p2, p1);
					clipping(_t, lightIntensity);
				}
				else if (t[outIndexes[0]].getY() > ScreenHeight() && t[outIndexes[1]].getY() > ScreenHeight())
				{
					Vec3 dirP1 = t[inIndex] - t[outIndexes[0]];
					dirP1 /= dirP1.getX();
					double mulP1 = (t[outIndexes[0]].getY() - ScreenHeight()) / dirP1.getY();
					p1 = t[outIndexes[0]] - (dirP1 * mulP1);
					p1.setY(ScreenHeight());

					Vec3 dirP2 = t[inIndex] - t[outIndexes[1]];
					dirP2 /= dirP2.getX();
					double mulP2 = (t[outIndexes[1]].getY() - ScreenHeight()) / dirP2.getY();
					p2 = t[outIndexes[1]] - (dirP2 * mulP2);
					p2.setY(ScreenHeight());

					Triangle _t = Triangle(t[inIndex], p2, p1);
					clipping(_t, lightIntensity);
				}
			}
			// If 1, delete one vertex, create two triangles
			// If 2, delete two vertices, create only one triangle
			// if 3, just delete the triangle
		}

		void drawProjection(Triangle& tProjection, double lightIntensity)
		{
		
			int outVertices = 0;
			for (unsigned int v = 0; v < 3; v++)
			{
				if (tProjection[v].getX() < 0.0 || tProjection[v].getX() > ScreenWidth() || tProjection[v].getY() < 0.0 || tProjection[v].getY() > ScreenHeight())
					outVertices++;
			}

			if (outVertices > 0)
			{
				clipping(tProjection , lightIntensity);
			}
			else
			{
				FillTriangle(
					tProjection[0].getX(),
					tProjection[0].getY(),
					tProjection[1].getX(),
					tProjection[1].getY(),
					tProjection[2].getX(),
					tProjection[2].getY(),
					olc::Pixel(255 * lightIntensity, 255 * lightIntensity, 255 * lightIntensity)
				);

				DrawTriangle(
						tProjection[0].getX(),
						tProjection[0].getY(),
						tProjection[1].getX(),
						tProjection[1].getY(),
						tProjection[2].getX(),
						tProjection[2].getY(),
						olc::GREEN
					);
			}
		}

		void generateProjection()
		{
			for (unsigned int i = 0; i < distanceBuffer.size(); i++)
			{
				Vec3 lightRotate = light;
				lightRotate.rotateX(-cameraRotation.getX());
				lightRotate.rotateY(-cameraRotation.getY());
				lightRotate.rotateZ(-cameraRotation.getZ());
				
				Vec3 lightDirection = distanceBuffer[i][0] - lightRotate + cameraPosition;
				lightDirection.normalize();

				double lightIntensity = -lightDirection.dot(distanceBuffer[i].normal());
				if (lightIntensity < 0.1)
					lightIntensity = 0.1;

				Triangle tProjection = Triangle(
					projection.mul(distanceBuffer[i][0]),
					projection.mul(distanceBuffer[i][1]),
					projection.mul(distanceBuffer[i][2]));

				for (unsigned int t = 0; t < 3; t++)
				{
					tProjection[t] = tProjection[t] + Vec3(1, 1, 0);

					tProjection[t].setX(tProjection[t].getX() * 0.5 * ScreenWidth());
					tProjection[t].setY(tProjection[t].getY() * 0.5 * ScreenHeight());
				}

				if (lightIntensity > 0)
				{
					drawProjection(tProjection, lightIntensity);
				}
			}
		}

	public : 
		
		GameEngine()
		{
			sAppName = "Easy Fast Game Engine";
			delta = 0.0;
			cameraPosition = Vec3(0, 0, -3);
			cameraRotation = Vec3(0, 0, 0);
			light = Vec3(-3, -10, 4);
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
			
			/*Mesh m = Mesh();
			
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

			Mesh prism = Mesh("simpleshape.obj");
			prism.translate(Vec3(0, 2, 0));

			addMesh(prism);*/

			Mesh test = Mesh();
			test.add(Triangle());
			test.translate(Vec3(0,0.5,0.5));
			test[0][1] = Vec3(-0.5,2,0);
			//test[0][0] = Vec3(0,0,0);
			//test[0][1] = Vec3(0, 0, 0);
			addMesh(test);

			/*Mesh teapot = Mesh("teapot.obj");
			teapot.translate(Vec3(0, 0, 5));
			addMesh(teapot);*/

			return true;
		}

		bool OnUserUpdate(float elapsedTime) override
		{
			FillRect(0, 0, ScreenWidth(), ScreenHeight(), olc::BLACK);

			delta += 1.0f * elapsedTime;

			// Camera Mouvement
			if (GetKey(olc::Key::UP).bPressed)
				cameraPosition += Vec3(0, 0, 1);

			if (GetKey(olc::Key::DOWN).bPressed)
				cameraPosition += Vec3(0, 0, -1);

			if (GetKey(olc::Key::LEFT).bPressed)
				cameraPosition += Vec3(-1, 0, 0);

			if (GetKey(olc::Key::RIGHT).bPressed)
				cameraPosition += Vec3(1, 0, 0);

			if (GetKey(olc::Key::SPACE).bPressed)
				cameraPosition += Vec3(0, -1, 0);

			if (GetKey(olc::Key::CTRL).bPressed)
				cameraPosition += Vec3(0, 1, 0);

			// Objects Mouvement
			if (GetKey(olc::Key::M).bPressed)
			{
				//for (unsigned int i = 0; i < mesh.size(); i++)
					mesh[0].translate(Vec3(1, 0, 0));
			}

			if (GetKey(olc::Key::K).bPressed)
			{
				//for (unsigned int i = 0; i < mesh.size(); i++)
					mesh[0].translate(Vec3(-1, 0, 0));
			}

			if (GetKey(olc::Key::O).bPressed)
			{
				//for (unsigned int i = 0; i < mesh.size(); i++)
				mesh[0].translate(Vec3(0, 0, 1));
			}

			if (GetKey(olc::Key::L).bPressed)
			{
				//for (unsigned int i = 0; i < mesh.size(); i++)
				mesh[0].translate(Vec3(0, 0, -1));
			}

			// Camera Rotation
			if (GetKey(olc::Key::Z).bPressed)
				cameraRotation -= Vec3(0.1, 0, 0);

			if (GetKey(olc::Key::S).bPressed)
				cameraRotation += Vec3(0.1, 0, 0);

			if (GetKey(olc::Key::Q).bPressed)
				cameraRotation += Vec3(0, 0.1, 0);

			if (GetKey(olc::Key::D).bPressed)
				cameraRotation -= Vec3(0, 0.1, 0);
			
			if (GetKey(olc::Key::A).bPressed)
				cameraRotation += Vec3(0, 0, 0.1);

			if (GetKey(olc::Key::E).bPressed)
				cameraRotation -= Vec3(0, 0, 0.1);

			for (unsigned int i = 0; i < mesh.size(); i++)
				show(mesh[i]);

			std::sort(distanceBuffer.begin(), distanceBuffer.end(), [](Triangle& t1, Triangle& t2)
			{
				Vec3 centerT1 = (t1[0] + t1[1] + t1[2]) / 3.0;
				Vec3 centerT2 = (t2[0] + t2[1] + t2[2]) / 3.0;

				double dist1 = centerT1.magnitude();
				double dist2 = centerT2.magnitude();

				return dist1 > dist2;
			});

			generateProjection();
			distanceBuffer.clear();

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
	if (ge.Construct(160, 160, 4, 4))
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