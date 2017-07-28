import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import math.sqrt
import math.cos
import math.sin

object MME{
  def main(args:Array[String]){
    println("Hello, World!")
    var arr1 = Nd4j.zeros(3)
    println("before",arr1)
    arr1 = matrot(arr1)
    println("after",arr1)
  }
  //Below are functions to change matrix dimensions
  //They use Mandel-Voigt notation to transform symmetric tensors
  def voigt_33_61(S33:INDArray): INDArray ={
    //Takes a matrix of dimensions (3 x 3 x 3 x 3) x N x M
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions (6 x 1 x 1 x 1) x N x M
    val transform_mat = Array(Array(0,1,2,1,2,0),
                              Array(0,1,2,2,0,1),
                              Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = S33.shape()
    var S61 = Nd4j.zeros(6,1,1,1,arr_size(4),arr_size(5))

    for(n<- 0 until (arr_size(4) - 1)){
      for(m <-0 until (arr_size(5) - 1)){
        for(i <- 0 until 5){
          val toPut = transform_mat(2)(i).asInstanceOf[Double]*S33.getDouble(transform_mat(0)(i).asInstanceOf[Int],transform_mat(1)(i).asInstanceOf[Int],0,0,n,m)
          S61.putScalar(Array(i,0,0,0,n,m),toPut)
        }
      }
    }
    S61
  }
  def voigt_333_36(d333:INDArray): INDArray ={
    //Takes a matrix of dimensions (3 x 3 x 3 x 1) x N x M
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions (3 x 6 x 1 x 1) x N x M
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = d333.shape()
    var d36 = Nd4j.zeros(3,6,1,1,arr_size(4),arr_size(5))

    for(n<- 0 until (arr_size(4) - 1)){
      for(m <-0 until (arr_size(5) - 1)){
        for(i <- 0 until 2){
          for(j<-0 until 5) {
            val toPut = transform_mat(2)(j).asInstanceOf[Double] * d333.getDouble(i,transform_mat(0)(j).asInstanceOf[Int], transform_mat(1)(j).asInstanceOf[Int], 0, n, m)
            d36.putScalar(Array(i, j, 0, 0, n, m), toPut)
          }
        }
      }
    }
    d36
  }
  def voigt_3333_66(C3333:INDArray): INDArray ={
    //Takes a matrix of dimensions (3 x 3 x 3 x 3) x N x M
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions (6 x 6 x 1 x 1) x N x M
    val transform_mat = Array(Array(0,1,2,1,2,0),
                              Array(0,1,2,2,0,1),
                              Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = C3333.shape()
    var C66 = Nd4j.zeros(6,6,1,1,arr_size(4),arr_size(5))

    for(n<- 0 until (arr_size(4) - 1)){
      for(m <-0 until (arr_size(5) - 1)){
        for(i<-0 until 5){
          for(j<-0 until 5) {
            val toPut = transform_mat(2)(i).asInstanceOf[Double] *transform_mat(2)(j).asInstanceOf[Double] * C3333.getDouble(transform_mat(0)(i).asInstanceOf[Int], transform_mat(1)(i).asInstanceOf[Int], transform_mat(0)(j).asInstanceOf[Int],transform_mat(1)(j).asInstanceOf[Int], n, m)
            C66.putScalar(Array(i, j, 0, 0, n, m), toPut)
          }
        }
      }
    }
    C66
  }
  def voigt_36_333(d36:INDArray): INDArray ={
    //takes a matrix of dimensions 3 x 6
    //and returns a matrix of dimensions 3 x 3 x 3
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    var d333 = Nd4j.zeros(3,3,3)

    for(i<- 0 until 2){
      for(j <-0 until 5) {
        val toPut = 1 / (transform_mat(2)(j).asInstanceOf[Double] * d36.getDouble(i, j))
        d333.putScalar(Array(i, transform_mat(0)(j).asInstanceOf[Int], transform_mat(1)(j).asInstanceOf[Int]), toPut)
        d333.putScalar(Array(i, transform_mat(1)(j).asInstanceOf[Int], transform_mat(0)(j).asInstanceOf[Int]), toPut)
      }
    }
    d333
  }
  def voigt_66_3333(C66:INDArray): INDArray ={
    //takes a matrix of dimensions 6 x 6
    //and returns a matrix of dimensions 3 x 3 x 3 x 3
    val transform_mat = Array(Array(0,1,2,1,2,0),
                              Array(0,1,2,2,0,1),
                              Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    var C3333 = Nd4j.zeros(3,3,3)

    for(t<- 0 until 5){
      for(u <-0 until 5) {
        val toPut = 1 / (transform_mat(2)(t).asInstanceOf[Double]*transform_mat(2)(u).asInstanceOf[Double]*C66.getDouble(t, u))
        C3333.putScalar(Array(transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(u).asInstanceOf[Int],transform_mat(1)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(t).asInstanceOf[Int], transform_mat(0)(u).asInstanceOf[Int],transform_mat(1)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(t).asInstanceOf[Int], transform_mat(1)(u).asInstanceOf[Int],transform_mat(0)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(u).asInstanceOf[Int],transform_mat(0)(u).asInstanceOf[Int]), toPut)
      }
    }
    C3333
  }
  //Below are rotations of matrices
  def rotate_ordre3(d333:INDArray,R:INDArray):INDArray ={
    //Takes a matrix of size (3 x 3 x 3 x 1) x N x M
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = d333.shape()
    var postRot = Nd4j.zeros(3,3,3,1,arr_size(4),arr_size(5))

    for(grain<-0 until (arr_size(4) - 1)){
      for(dir<-0 until (arr_size(5) - 1)){
        for(i<-0 until 2){
          for(j<-0 until 2){
            for(k<-0 until 2){
              for(m<-0 until 2){
                for(n<-0 until 2){
                  for(o<- 0 until 2){
                    var toPut = postRot.getDouble(i,j,k,0,grain,dir)+R.getDouble(i,m,0,0,grain,dir)*R.getDouble(j,n,0,0,grain,dir)*R.getDouble(k,o,0,0,grain,dir)*d333.getDouble(m,n,o,0,grain,dir)
                    postRot.putScalar(Array(i,j,k,0,grain,dir),toPut)
                  }
                }
              }
            }
          }
        }
      }
    }
    postRot
  }
  def rotate_ordre4(C3333:INDArray,R:INDArray):INDArray ={
    //Takes a matrix of size (3 x 3 x 3 x 3) x N x M
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = C3333.shape()
    var postRot = Nd4j.zeros(3,3,3,3,arr_size(4),arr_size(5))

    for(grain<-0 until (arr_size(4) - 1)){
      for(dir<-0 until (arr_size(5) - 1)){
        for(i<-0 until 2){
          for(j<-0 until 2){
            for(k<-0 until 2){
              for(l<-0 until 2) {
                for (m <-0 until 2) {
                  for (n <-0 until 2) {
                    for (o <-0 until 2) {
                      for (p <- 0 until 2) {
                        var toPut = postRot.getDouble(i, j, k, l, grain, dir) + R.getDouble(i, m, 0, 0, grain, dir) * R.getDouble(j, n, 0, 0, grain, dir) * R.getDouble(k, o, 0, 0, grain, dir)*R.getDouble(l,p,0,0,grain,dir) * C3333.getDouble(m, n, o, p, grain, dir)
                        postRot.putScalar(Array(i, j, k, 0, grain, dir), toPut)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    postRot
  }
  def matrot(angles:INDArray): INDArray ={
    //converts a 3 x 1 orientation vector into a 3 x 3 matrix
    //UNSURE HOW DIMENSIONALITY OF ANGLES ARRAY WILL LINE UP IN DONN_POLYC BE CAREFUL
    var A_rot = Nd4j.zeros(3,3)

    val c1 = cos(angles.getDouble(0))
    val s1 = sin(angles.getDouble(0))
    val c2 = cos(angles.getDouble(1))
    val s2 = sin(angles.getDouble(1))
    val c3 = cos(angles.getDouble(2))
    val s3 = sin(angles.getDouble(2))

    A_rot.putScalar(Array(0,0) ,c1*c3-s1*c2*s3)
    A_rot.putScalar(Array(0,1),-c1*s3-c2*c3*s1)
    A_rot.putScalar(Array(0,2),s1*s2)
    A_rot.putScalar(Array(1,0),c3*s1+c1*c2*s3)
    A_rot.putScalar(Array(1,1),c1*c2*c3-s1*s3)
    A_rot.putScalar(Array(1,2),-c1*s2)
    A_rot.putScalar(Array(2,0),s2*s3)
    A_rot.putScalar(Array(2,1),c3*s2)
    A_rot.putScalar(Array(2,2),c2)

    A_rot
  }
  //below are the eshelby funtions
  def Eshelby33_muphylin(C_inf:INDArray,form_inclu:Int):INDArray = {
    var Ndemag33:INDArray = Nd4j.zeros(3)
    if(form_inclu == 1){
      Ndemag33.muli(1/3)
    }else{
      Ndemag33.putScalar(0,0,1/2)
      Ndemag33.putScalar(1,1,1/2)
    }
    Ndemag33
  }
  def Eshelby33_muphylin(C_inf:INDArray,icalc:Int,semi_axis:Array[Array[Int]]):INDArray = {
    var Ndemag33:INDArray = Nd4j.zeros(3)
    if(icalc == 1){
      Ndemag33.muli(1/3)
    }else{
      Ndemag33.putScalar(0,0,1/2)
      Ndemag33.putScalar(1,1,1/2)
    }
    Ndemag33
  }
}
