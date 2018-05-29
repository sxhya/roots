import org.junit.Test

import java.util.Arrays

/**
 * Created by user on 6/15/17.
 */

class Tests {
    @Test
    fun test() {
        val k = arrayOf(longArrayOf(4, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0))

        val im = IntegerMatrix(k)
        val computer = IntegerMatrix.SNFComputer(im)
        computer.doComputeSNF()
        println(Arrays.toString(computer.s.getDiagonalEntries()))
    }

    @Test
    fun testInverse() {
        val m = Matrix.random(4, 5.0, false)
        val mi = m.inverse()
        println(m)
        println()
        println(mi)
        println()
        println(m.mul(mi))
    }

}
