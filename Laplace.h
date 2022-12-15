#include <iostream>
#include <vector>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "Eigen/Core"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
struct simpleTraits : public OpenMesh::DefaultTraits {
    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color);
    EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::TriMesh_ArrayKernelT<simpleTraits> simpleMesh;

#ifndef OPENMESHPOLYSCOPE_LAPLACE_H
#define OPENMESHPOLYSCOPE_LAPLACE_H

////Fonction permettant de faire i iteration de lissage laplacien du maillage mesh
simpleMesh Lissage(int iteration, simpleMesh &mesh);

//// n : taille du n-ring
//// indice_point : numéro de l'indice du vertex centrale du ring
//// mesh : le maillage sur lequel appliqués la fonction
//// value : correspond à ma custom properties
std::vector<MyMesh::VertexHandle> n_ring(int n, int indice_point, simpleMesh &mesh,OpenMesh::VPropHandleT<double> value);

////Calcule le 1-ring d'un vertex
//// mesh : le maillage sur lequel appliqués la fonction
//// vh : le centre sur lequel on va calculer le 1-ring
std::vector<MyMesh::VertexHandle> premier_ring(MyMesh::VertexHandle vh,simpleMesh &mesh);

////Résout le problème AX = B et remplis notre custom property avec les valeurs de X
//// ring : listes des vertexHandle des points appartenants a un ring
//// mesh : le maillage sur lequel appliqués la fonction
//// value : correspond à ma custom properties
void create_Eigen_values(std::vector<MyMesh::VertexHandle> ring,simpleMesh &mesh,OpenMesh::VPropHandleT<double> value);

#endif //OPENMESHPOLYSCOPE_LAPLACE_H
