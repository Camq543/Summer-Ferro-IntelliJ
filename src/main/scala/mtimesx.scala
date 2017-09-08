import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j

object mtimesx {
  //Takes two 4-dimensional arrays of size i x j x k x l and i x j x l x n, and a char denoting the mode
  //returns an array of size i x j x k x n which is found by multiplying all combinations of two-D arrays of size k x l and l x n
  def mul(arr1:INDArray,arr2:INDArray,mode:Char): INDArray ={
    var outArr :INDArray = null
    val maxGrain = math.max(arr1.shape()(0),arr2.shape()(0))
    val maxDir = math.max(arr1.shape()(1),arr2.shape()(1))
    var grain1 = 0
    var dir1 = 0
    var grain2 = 0
    var dir2 = 0
    if( mode == 'N') {
      if (arr1.shape()(3) != arr2.shape()(2)) {
        println("ERROR in mtimesx, row mismatch")
      }
      outArr = Nd4j.create(maxGrain,maxDir,arr1.shape()(2),arr2.shape()(3))
    }else if (mode =='T') {
      if (arr1.shape()(3) != arr2.shape()(3)) {
        println("ERROR in mtimesx, row mismatch")
      }
      outArr = Nd4j.create(maxGrain,maxDir,arr1.shape()(2),arr2.shape()(2))
    }

    var toPut:INDArray = Nd4j.create(arr1.shape()(2),arr2.shape()(3))
    for(grain <- 0 until maxGrain){
      for(dir <- 0 until maxDir){
        grain1 = grain
        grain2 = grain
        dir1 = dir
        dir2 = dir
        if(grain>=arr1.shape()(0)){
          grain1 = arr1.shape()(0) - 1
        }
        if(dir>=arr1.shape()(1)){
          dir1 = arr1.shape()(1) - 1
        }
        if(grain>=arr2.shape()(0)){
          grain2 = arr2.shape()(0) - 1
        }
        if(dir>=arr2.shape()(1)){
          dir2 = arr2.shape()(1) - 1
        }
        if(mode == 'N') {
          //println("grain",grain,"dir",dir,"grain1",grain1,"dir1",dir1,"grain2",grain2,"dir2",dir2)
          //println("arr1",arr1.getRow(grain1).getRow(dir1),"arr2",arr2.getRow(grain2).getRow(dir2))
          toPut = arr1.getRow(grain1).getRow(dir1).mmul(arr2.getRow(grain2).getRow(dir2))
        } else if(mode =='T') {
          toPut = arr1.getRow(grain1).getRow(dir1).mmul(arr2.getRow(grain2).getRow(dir2).transpose())
        }
        outArr.getRow(grain).putRow(dir,toPut)
      }
    }
    outArr
  }
}
