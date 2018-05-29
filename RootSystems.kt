/**
 * Created by user on 1/17/15.
 */

class RootSystems {
    companion object {

        val e8base = Matrix(arrayOf(doubleArrayOf(-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5),
                doubleArrayOf( 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0)))

        val e7base = Matrix(arrayOf(doubleArrayOf(-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, Math.sqrt(2.0) / 2),
                doubleArrayOf( 0.0,  0.0,  0.0,  0.0,  1.0, -1.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  0.0,  1.0,  1.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  1.0, -1.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  1.0, -1.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  1.0, -1.0,  0.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 1.0, -1.0,  0.0,  0.0,  0.0,  0.0, 0.0)))

        val e6base = Matrix(arrayOf(doubleArrayOf( 1.0, -1.0,  0.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  1.0, -1.0, 0.0),
                doubleArrayOf( 0.0,  1.0, -1.0,  0.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  1.0, -1.0,  0.0, 0.0),
                doubleArrayOf( 0.0,  0.0,  0.0,  1.0,  1.0, 0.0),
                doubleArrayOf(-0.5, -0.5, -0.5, -0.5, -0.5, Math.sqrt(3.0) / 2)))

        val f4base = Matrix(arrayOf(doubleArrayOf( 1.0, -1.0,  0.0,  0.0),
                doubleArrayOf( 0.0,  1.0, -1.0,  0.0),
                doubleArrayOf( 0.0,  0.0,  1.0,  0.0),
                doubleArrayOf(-0.5, -0.5, -0.5, -0.5)))


        fun dlbase(l: Int): Matrix {
            val coos = Array(l) { DoubleArray(l) }
            for (i in 0 until l) {
                coos[i][i] = 1.0
                if (i < l - 1) coos[i][i + 1] = -1.0
            }
            coos[l - 1][l - 2] = 1.0
            return Matrix(coos)
        }

        fun albase(l: Int): Matrix {
            val coos = Array(l) { DoubleArray(l) }
            for (i in 0 until l - 1) {
                coos[i][i] = 1.0
                coos[i][i + 1] = -1.0
            }
            val b = (l - 1 + Math.sqrt((l + 1).toDouble())) / l
            for (i in 0 until l - 1) coos[l - 1][i] = b - 1
            coos[l - 1][l - 1] = b
            return Matrix(coos)
        }

        fun blbase(l: Int): Matrix {
            val coos = clbase(l).myCoo
            coos[l - 1][l - 1] = 1.0
            return Matrix(coos)
        }


        fun clbase(l: Int): Matrix {
            val coos = Array(l) { DoubleArray(l) }
            for (i in 0 until l - 1) {
                coos[i][i] = 1.0
                coos[i][i + 1] = -1.0
            }
            for (i in 0 until l - 2) coos[l - 1][i] = 0.0
            coos[l - 1][l - 1] = 2.0
            return Matrix(coos)
        }

    }
}

