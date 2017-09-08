import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import org.nd4j.linalg.inverse

class polyc(monoCrystal:monoc) {
  //sets up the polycrystal
  //n1, n2, and n3 are used to create a distribution of orientations
  //they correspond to the x, y, and z dimensions
  val n1:Int = 13
  val n2:Int = 7
  val n3:Int = 6
  val nbg = n1*n2*n3
  var angles:INDArray = Nd4j.create(nbg,1,3,1)
  anglesetter()
  //angles array is set up, now we turn it into an array of grain orientation matrixes
  var M_rot:INDArray = Nd4j.create(nbg,1,3,3)
  for(i<-0 until nbg){
    var toPut = MME.matrot(angles.getRow(i).getRow(0))
    M_rot.getRow(i).putRow(0,toPut)
  }

  //transfer over the monocrystal data, and put the grains into alignment with monocrystal orientation
  var dir100:INDArray = mtimesx.mul(M_rot,monoCrystal.dir100,'N')
  var dir110:INDArray = mtimesx.mul(M_rot,monoCrystal.dir110,'N')
  var dir111:INDArray = mtimesx.mul(M_rot,monoCrystal.dir111,'N')

  var Pol:INDArray = dir100.mul(monoCrystal.Pol0)
  var epsferro61:INDArray = MME.voigt_33_61(mtimesx.mul(mtimesx.mul(M_rot,monoCrystal.epsferro33,'N'),M_rot,'T'))
  var kappa33:INDArray = mtimesx.mul(mtimesx.mul(M_rot,monoCrystal.kappa33,'N'),M_rot,'T')
  var d36:INDArray = MME.voigt_333_36(dHelper())
  var C66:INDArray = MME.voigt_3333_66(cHelper())
  var invC66:INDArray = Nd4j.create(nbg,monoCrystal.nbdom,6,6)
  for(n<-0 until nbg){
    for(m<-0 until monoCrystal.nbdom){
      invC66.getRow(n).putRow(m,inverse.InvertMatrix.invert(C66.getRow(n).getRow(m),false))
    }
  }


  println("Polycrystal setup finished")

  def anglespace(min:Double, max:Double, len:Int): INDArray ={
    //uses linspace to give an array with evenly distributed values across a range
    //cuts out last element
    var out = Nd4j.create(len)
    var space = Nd4j.linspace(min, max, len + 1)
    for(i <- 0 until len){
      out.putScalar(i,space.getDouble(i))
    }
    out
  }
  def anglesetter(): Unit ={
    //constructs the array which holds the orientation of all polycrystal grains
    val angle1 = anglespace(0,2*math.Pi,n1)
    val angle2 = anglespace(0,1,n2)
    for(i<-0 until n2){
      angle2.putScalar(i,math.acos(angle2.getDouble(i)))
    }
    val angle3 = anglespace(0,math.Pi/2,n3)
    var count = 0
    var toPut = Nd4j.create(3)
    for(i<- 0 until n1){
      for(j<- 0 until n2){
        for(k<- 0 until n3){
          toPut = Nd4j.create(Array(angle1.getDouble(i),angle2.getDouble(j),angle3.getDouble(k)))
          this.angles.getRow(count).getRow(0).putColumn(0,toPut)
          count += 1
        }
      }
    }
  }
  def dHelper():INDArray = {
    var out: INDArray = null
    var dFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3,3)
    var rotFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3)
    for(i<-0 until nbg){
      dFill.putRow(i,monoCrystal.dmonoc333.getRow(0))
      for(j<-0 until monoCrystal.nbdom){
        rotFill.getRow(i).putRow(j,M_rot.getRow(i).getRow(0))
      }
    }
    out = MME.rotate_ordre3(dFill,rotFill)
    out
  }
  def cHelper(): INDArray ={
    var out: INDArray = null
    var cFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3,3,3)
    var rotFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3)
    for(i<-0 until nbg){
      cFill.putRow(i,monoCrystal.Cmonoc3333.getRow(0))
      println(cFill.shape()(4))
      for(j<-0 until monoCrystal.nbdom){
        rotFill.getRow(i).putRow(j,M_rot.getRow(i).getRow(0))
      }
    }
    println("cFill",cFill)
    out = MME.rotate_ordre4(cFill,rotFill)
    println(out)
    out
  }
}
