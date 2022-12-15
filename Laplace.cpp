#include <iostream>
#include <vector>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
using namespace std;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
struct simpleTraits : public OpenMesh::DefaultTraits {
    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color);
    EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::TriMesh_ArrayKernelT<simpleTraits> simpleMesh;

using Eigen::MatrixXd;



simpleMesh Lissage(int iteration, simpleMesh &mesh)
{

    // this vector stores the computed centers of gravity
    std::vector<MyMesh::Point>  cogs;
    std::vector<MyMesh::Point>::iterator cog_it;
    cogs.reserve(mesh.n_vertices());
    // smoothing mesh argv[1] times
    MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());
    MyMesh::VertexVertexIter    vv_it;
    MyMesh::Point               cog;
    MyMesh::Scalar              valence;
    unsigned int                i;
    float alpha = 0.9;
    for (i=0; i < iteration; i++)
    {
        cogs.clear();
        for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
        {
            cog[0] = cog[1] = cog[2] = valence = 0.0;

            for (vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it)
            {
                cog += mesh.point( *vv_it );
                ++valence;
            }
            MyMesh::Point v_i = mesh.point(*v_it);
            cogs.push_back(v_i*alpha+(1-alpha)*(cog / valence));
        }

        for (v_it=mesh.vertices_begin(), cog_it=cogs.begin();
             v_it!=v_end; ++v_it, ++cog_it)
            if ( !mesh.is_boundary( *v_it ) )
                mesh.set_point( *v_it, *cog_it );
    }

    return mesh;
}

std::vector<MyMesh::VertexHandle> premier_ring(MyMesh::VertexHandle vh,simpleMesh &mesh){

    std::vector<MyMesh::VertexHandle> ring;
    ring.clear();

    MyMesh::VertexVertexIter    vv_it;
    for (vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
        ring.push_back(*vv_it);
    }
    return ring;

}


std::vector<MyMesh::VertexHandle> n_ring(int n, int indice_point, simpleMesh &mesh,OpenMesh::VPropHandleT<double> value){

////    mes anneaux
    std::vector<MyMesh::VertexHandle>  ring;
    std::vector<MyMesh::VertexHandle> anno;
    std::vector<MyMesh::VertexHandle> annok;
    std::vector<MyMesh::VertexHandle> anno_tmp;


    auto N = static_cast<double>(n);


    MyMesh::VertexHandle point;
    MyMesh::VertexHandle centre = mesh.vertex_handle(indice_point);
    anno.push_back(centre);

    for(int i = 0 ; i<=n ; i++){ // i == iter sur le nombre de ring
        for(auto j : anno) { // pour visiter les elem de anno
            anno_tmp.clear();
            point = j;
            if (mesh.property(value, point) == 0.0) {
                auto I = static_cast<double>(i);
                mesh.property(value, point) = double (1.0*(N - I) / N);
//                cout << mesh.property(value, point) << "\n";
                anno_tmp.push_back(point);
                ring.push_back(point);
            }
        }
        anno.clear();
        for(auto k : anno_tmp){
            point = k;
            annok = premier_ring(point,mesh);
            for(auto l : annok){
                anno.push_back(l);
            }
        }
    }

    return ring;
}




void create_Eigen_values(std::vector<MyMesh::VertexHandle> ring,simpleMesh &mesh,OpenMesh::VPropHandleT<double> value){

    ////On réccupère la taille de notre vecteur et on initialise nos matrices/vecteurs pour les remplir
    int nb_v = ring.size();
    Eigen::MatrixXd A(nb_v,nb_v);
    Eigen::VectorXd B(nb_v);
    MyMesh::VertexHandle v_act;
    double val,valence;
    std::vector<int> index_v;


    for(int i = 0 ; i<nb_v ; i++){
        for(int j = 0 ; j<nb_v ; j++){
            A(i,j) = 0.;
        }
        B(i) = 0.;
    }

    for(int i = 0 ; i<nb_v ; i++){
        v_act = ring.at(i);
        val = mesh.property(value, v_act);
        if(val == 1.) {
            B(i) = 1.;
            A(i, i) = 1.0;
        }
        else
            {
                if (val == 0.) {
                    A(i, i) = 1.;
                } else {
                    A(i, i) = -1.0;


                    vector<MyMesh::VertexHandle> anno1 = premier_ring(v_act, mesh);
                    for (auto vh: anno1) {
                        for (int k = 0; k < ring.size(); k++) {
                            valence = anno1.size();
                            if (vh == ring.at(k)) {
                                A(i, k) = 1.0 / valence;
                            }
                        }
                    }
                }
            }

    }
    Eigen::VectorXd X;
    X = A.colPivHouseholderQr().solve(B);
//    cout << X << endl;
    for (int i = 0 ; i<ring.size() ; i++ ){
        v_act = ring.at(i);
        mesh.property(value,v_act) = X(i);
    }

}
