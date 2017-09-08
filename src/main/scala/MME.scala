import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import math.sqrt
import math.cos
import math.sin

object MME{

  def main(args:Array[String]){
    var monoCrystal = new monoc()
    var polyCrystal = new polyc(monoCrystal)
    var arr1 = Nd4j.create(6,6)
    var arr2 = arr1.repmat(3,3,6,6)
  }
  //Set-up functions
  //Monocrystal set-up
  def donn_monoc_tetra(): Unit ={

  }
  //Below are functions to change matrix dimensions
  //They use Mandel-Voigt notation to transform symmetric tensors
  def voigt_33_61(S33:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (6 x 1)
    val transform_mat = Array(Array(0,1,2,1,2,0),
                              Array(0,1,2,2,0,1),
                              Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = S33.shape()
    var S61 = Nd4j.zeros(arr_size(0),arr_size(1),6,1)


    for(n<- 0 until arr_size(0)){
      for(m <-0 until arr_size(1)){
        for(i <- 0 until 6){
          val toPut = transform_mat(2)(i).asInstanceOf[Double]*S33.getDouble(n,m,transform_mat(0)(i).asInstanceOf[Int],transform_mat(1)(i).asInstanceOf[Int])
          S61.putScalar(Array(n,m,i,0),toPut)
        }
      }
    }
    S61
  }
  def voigt_333_36(d333:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (3 x 6)
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = d333.shape()
    var d36 = Nd4j.zeros(arr_size(0),arr_size(1),3,6)

    for(n<- 0 until (arr_size(0))){
      for(m <-0 until (arr_size(1))){
        for(i <- 0 until 3){
          for(j<-0 until 6) {
            val toPut = transform_mat(2)(j).asInstanceOf[Double] * d333.getDouble(n,m,i,transform_mat(0)(j).asInstanceOf[Int], transform_mat(1)(j).asInstanceOf[Int])
            d36.putScalar(Array(n,m,i,j), toPut)
          }
        }
      }
    }
    d36
  }
  def voigt_3333_66(C3333:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (6 x 6)
    val transform_mat = Array(Array(0,1,2,1,2,0),
                              Array(0,1,2,2,0,1),
                              Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = C3333.shape()
    var C66 = Nd4j.zeros(arr_size(0),arr_size(1),6,6)

    for(n<- 0 until arr_size(0)){
      for(m <-0 until arr_size(1)){
        for(i<-0 until 6){
          for(j<-0 until 6) {
            val toPut = transform_mat(2)(i).asInstanceOf[Double] *transform_mat(2)(j).asInstanceOf[Double] * C3333.getDouble(n,m,transform_mat(0)(i).asInstanceOf[Int], transform_mat(1)(i).asInstanceOf[Int], transform_mat(0)(j).asInstanceOf[Int],transform_mat(1)(j).asInstanceOf[Int])
            C66.putScalar(Array(n,m,i,j), toPut)
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

    for(i<- 0 until 3){
      for(j <-0 until 6) {
        val toPut = (1 / (transform_mat(2)(j).asInstanceOf[Double]) * d36.getDouble(i, j))
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
    var C3333 = Nd4j.zeros(3,3,3,3)

    for(t<- 0 until 6){
      for(u <-0 until 6) {
        val toPut = (1 / (transform_mat(2)(t).asInstanceOf[Double]*transform_mat(2)(u).asInstanceOf[Double])*C66.getDouble(t, u))
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
    //Takes a matrix of size N x M x (3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = d333.shape()
    var postRot = Nd4j.zeros(arr_size(0),arr_size(1),3,3,3)

    for(grain<-0 until arr_size(0)){
      for(dir<-0 until arr_size(1)){
        for(i<-0 until 3){
          for(j<-0 until 3){
            for(k<-0 until 3){
              for(m<-0 until 3){
                for(n<-0 until 3){
                  for(o<- 0 until 3){
                    var toPut = postRot.getDouble(grain,dir,i,j,k)+R.getDouble(grain,dir,i,m)*R.getDouble(grain,dir,j,n)*R.getDouble(grain,dir,k,o)*d333.getDouble(grain,dir,m,n,o)
                    postRot.putScalar(Array(grain,dir,i,j,k),toPut)
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
    //Takes a matrix of size N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = C3333.shape()
    var postRot = Nd4j.zeros(arr_size(0),arr_size(1),3,3,3,3)

    for(grain<-0 until arr_size(0)){
      for(dir<-0 until arr_size(1)){
        for(i<-0 until 3){
          for(j<-0 until 3){
            for(k<-0 until 3){
              for(l<-0 until 3) {
                for (m <-0 until 3) {
                  for (n <-0 until 3) {
                    for (o <-0 until 3) {
                      for (p <- 0 until 3) {
                        var toPut = postRot.getDouble(grain, dir, i, j, k, l) + (R.getDouble(grain, dir, i, m) * R.getDouble(grain, dir, j, n) * R.getDouble(grain, dir, k, o)*R.getDouble(grain, dir, l, p) * C3333.getDouble(grain, dir, m, n, o, p))
                        postRot.putScalar(Array(grain, dir, i, j, k, l), toPut)
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
  //below are the eshelby functions
  def Eshelby33_muphylin(C_inf:INDArray,form_inclu:Int):INDArray = {
    var Ndemag33:INDArray = Nd4j.zeros(3,3)
    if(form_inclu == 1){
      Ndemag33 = Nd4j.eye(3)
      Ndemag33.muli(1.0/3)
    }else{
      Ndemag33.putScalar(0,0,1.0/2)
      Ndemag33.putScalar(1,1,1.0/2)
    }
    Ndemag33
  }
  def Eshelby66_muphylin(C_inf:INDArray,icalc:Int,semi_axis:Array[Array[Int]]):INDArray = {
    var SEsh66:INDArray = Nd4j.zeros(3)
    if (icalc != 1){
      println("ERROR IN ESHELBY66: ICALC INPUT WRONG")
    }else{
      var nu0 = C_inf.getDouble(0,1)/(C_inf.getDouble(0,0) + C_inf.getDouble(0,1))
      val ones = Nd4j.ones(3,3)
      val zeros = Nd4j.zeros(3,3)
      val eye = Nd4j.eye(3)
      var JJ = Nd4j.concat(0,Nd4j.concat(1,ones,zeros),Nd4j.concat(1,zeros,zeros)).mul(1.0/3)
      var KK = Nd4j.concat(0,Nd4j.concat(1,ones.mul(-1).add(eye.mul(3)),zeros),Nd4j.concat(1,zeros,eye.mul(3))).mul(1.0/3)
      SEsh66 = JJ.mul((1+nu0)/(3.0*(1-nu0))).add(KK.mul((2/15)*(4-5*nu0)/(1-nu0)))
    }
    SEsh66
  }
}
