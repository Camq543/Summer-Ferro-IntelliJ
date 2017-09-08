//import org.nd4j.linalg.api.complex.IComplexNDArray
import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j

//class to hold monocrystal
//performs setup from donn_monoc_tetra in constructor
class monoc {

  //constants
  val nbdom = 6
  val AS:Double = "1E-5".toDouble
  val Pol0: Double = 0.3
  //set the monocrystal direction matrixes
  //these will be used to set up the polycrystal grains
  var dir100:INDArray = Nd4j.zeros(1,6,3,1)
  dir100.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(0.0),Array(0.0))))
  dir100.getRow(0).putRow(1,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(0.0))))
  dir100.getRow(0).putRow(2,Nd4j.create(Array(Array(0.0),Array(1.0),Array(0.0))))
  dir100.getRow(0).putRow(3,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(0.0))))
  dir100.getRow(0).putRow(4,Nd4j.create(Array(Array(0.0),Array(0.0),Array(1.0))))
  dir100.getRow(0).putRow(5,Nd4j.create(Array(Array(0.0),Array(0.0),Array(-1.0))))
  var dir110:INDArray = Nd4j.zeros(1,12,3,1)
  dir110.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(1.0),Array(0.0))))
  dir110.getRow(0).putRow(1,Nd4j.create(Array(Array(1.0),Array(-0.0),Array(0.0))))
  dir110.getRow(0).putRow(2,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(0.0))))
  dir110.getRow(0).putRow(3,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(0.0))))
  dir110.getRow(0).putRow(4,Nd4j.create(Array(Array(0.0),Array(1.0),Array(1.0))))
  dir110.getRow(0).putRow(5,Nd4j.create(Array(Array(0.0),Array(1.0),Array(-1.0))))
  dir110.getRow(0).putRow(6,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(1.0))))
  dir110.getRow(0).putRow(7,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(-1.0))))
  dir110.getRow(0).putRow(8,Nd4j.create(Array(Array(1.0),Array(0.0),Array(1.0))))
  dir110.getRow(0).putRow(9,Nd4j.create(Array(Array(1.0),Array(0.0),Array(-1.0))))
  dir110.getRow(0).putRow(10,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(1.0))))
  dir110.getRow(0).putRow(11,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(-1.0))))
  var dir111:INDArray = Nd4j.ones(1,8,3,1)
  dir111.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(1.0),Array(1.0))))
  dir111.getRow(0).putRow(1,Nd4j.create(Array(Array(1.0),Array(1.0),Array(-1.0))))
  dir111.getRow(0).putRow(2,Nd4j.create(Array(Array(1.0),Array(-1.0),Array(1.0))))
  dir111.getRow(0).putRow(3,Nd4j.create(Array(Array(1.0),Array(-1.0),Array(-1.0))))
  dir111.getRow(0).putRow(4,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(1.0))))
  dir111.getRow(0).putRow(5,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(-1.0))))
  dir111.getRow(0).putRow(6,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(1.0))))
  dir111.getRow(0).putRow(7,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(-1.0))))
  dir111.muli(1.0/math.sqrt(1/3))
  //other arrays
  var dir_ref: INDArray = Nd4j.create(Array(0,0,1.0)).reshape(3,1)
  var Pol: INDArray = this.dir100.mul(this.Pol0)
  def epsferroHelper():INDArray = {
    //sets up the epsferro33 array
    val epz0 = "2E-3".toDouble
    val mat_id = Nd4j.create(1,nbdom,3,3)
    val temp = Nd4j.eye(3)
    for(i <- 0 until nbdom) {
      mat_id.getRow(0).putRow(i, temp)
    }
    mtimesx.mul(this.dir100,this.dir100,'T').mul(3).sub(mat_id).mul(epz0/2)
  }
  var epsferro33: INDArray = epsferroHelper()

  //setup Tensors
  var Cmonoc3333: INDArray = Nd4j.zeros(1,nbdom,3,3,3,3)
  var kappa33: INDArray = Nd4j.zeros(1,nbdom,3,3)
  var dmonoc333: INDArray = Nd4j.zeros(1,nbdom,3,3,3)

  def tensorHelper(): INDArray ={
    //finds the Rotation matrix R which is used to calculate cmonoc3333, kappa33, and dmonoc333


    //math below is to calculate the "R" matrix
    "ISSUES"
    var R = Nd4j.zeros(1,nbdom,3,3)
    var cos_t = 1.0
    var dir_r = Nd4j.create(Array(0.0,0,0)).reshape(3,1)
    var mat_cross = Nd4j.zeros(3,3)
    var toPut = Nd4j.create(3,3)

    for(i <- 0 until nbdom){
      cos_t = dir_ref.transpose().mmul(this.dir100.getRow(0).getRow(i)).getDouble(0)
      dir_r.putColumn(0,Nd4j.create(Array(this.dir100.getDouble(0,i,1,0),-1*this.dir100.getDouble(0,i,0,0),0)))
      if(dir_r.norm2Number().asInstanceOf[Double] != 0){
        dir_r.divi(dir_r.norm2Number())
      }

      mat_cross.putRow(0,Nd4j.create(Array(0,0,dir_r.getDouble(1,0))))
      mat_cross.putRow(1,Nd4j.create(Array(0,0,-1*dir_r.getDouble(0,0))))
      mat_cross.putRow(2,Nd4j.create(Array(-dir_r.getDouble(1,0),dir_r.getDouble(0,0),0)))

    

      toPut = Nd4j.eye(3).mul(cos_t).add(dir_r.mmul(dir_r.transpose()).mul(1-cos_t)).sub(mat_cross.mul(math.sqrt(1-(cos_t * cos_t))))
      R.getRow(0).putRow(i,toPut)
    }
    R
  }
  var R:INDArray =tensorHelper()


  //now, with R, we calculate Cmonoc, kappa, and dmonoc, which are used later
  //array setup
  //constants necessary for setup
  def CmonocHelper(): INDArray = {
    val C11 = "106E9".toDouble
    val C12 = "62e9".toDouble
    val C44 = "44E9".toDouble
    var Cmonoc = Nd4j.create(Array(Array(C11, C12, C12, 0, 0, 0),
      Array(C11, C11, C12, 0, 0, 0),
      Array(C12, C12, C11, 0, 0, 0),
      Array(0, 0, 0, C44, 0, 0),
      Array(0, 0, 0, 0, C44, 0),
      Array(0, 0, 0, 0, 0, C44)))
    Cmonoc
  }
  var Cmonoc = CmonocHelper()

  def dmonocHelper():INDArray = {
    val d31 = "-2.1E-10".toDouble
    val d33 = "4.5E-10".toDouble
    val d15 = "5.8E-10".toDouble
    var dmonoc = Nd4j.create(Array(Array(0, 0, 0, 0, d15, 0),
      Array(0, 0, 0, d15, 0, 0),
      Array(d31, d31, d33, 0, 0, 0)))
    dmonoc
  }
  var dmonoc = dmonocHelper()

  def kappaHelper(): INDArray ={
    val k33 = "2.00E-8".toDouble
    val k11 = "2.00E-8".toDouble
    var kappa = Nd4j.create(Array(Array(k11,0,0),
      Array(0,k11,0),
      Array(0,0,k33)))
    kappa
  }
  var kappa = kappaHelper()

  var Cmonoc3333temp:INDArray = MME.voigt_66_3333(Cmonoc)
  var dmonoc333temp:INDArray = MME.voigt_36_333(dmonoc)

  //with each of the tensors calculated, we create arrays of N x M x T where
  //T is the tensor matrix Cmonoc, kappa, or Dmonoc

  for(i<- 0 until nbdom){
    this.Cmonoc3333.getRow(0).putRow(i,Cmonoc3333temp)
  }
  this.Cmonoc3333 = MME.rotate_ordre4(this.Cmonoc3333,R)
  for(i<-0 until nbdom){
    this.kappa33.getRow(0).putRow(i,kappa)
  }
  this.kappa33 = mtimesx.mul(mtimesx.mul(this.kappa33,R,'N'),R,'T')
  for(i<-0 until nbdom){
    this.dmonoc333.getRow(0).putRow(i,dmonoc333temp)
  }
  this.dmonoc333 = MME.rotate_ordre3(this.dmonoc333,R)

  println("Monocrystal setup finished")

}
