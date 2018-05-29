import java.util.ArrayList
import java.util.Arrays

class GCMExperiments {
    companion object {
        private fun convertGCMtoAwt(gcm: IntegerMatrix): IntegerMatrix {
            val n = gcm.height()
            assert(n == gcm.width())
            val result = IntegerMatrix(Array((n * n - n) / 2) { LongArray(n) })
            var counter = 0
            for (i in 0 until n)
                for (j in i + 1 until n) {
                    result.myCoo[counter][i] = gcm.myCoo[j][i]
                    result.myCoo[counter][j] = -gcm.myCoo[i][j]
                    counter++
                }
            return result
        }

        private fun gcmexperiment(name: String, gcm: Array<LongArray>) {
            val oddColumns = ArrayList<Int>()
            val evenColumns = ArrayList<Int>()

            val cm = IntegerMatrix(gcm)
            assert(cm.width() == cm.height())

            for (i in 0 until cm.width()) {
                var even = true
                for (j in 0 until cm.height())
                    if (cm.myCoo[j][i] % 2 != 0L)
                        even = false
                (if (even) evenColumns else oddColumns).add(i)
            }

            val p = IntArray(cm.width())
            var c = 0
            for (oc in oddColumns) p[c++] = oc
            for (oc in evenColumns) p[c++] = oc
            // o, e, o --> p = 0, 2 | 1


            val gcm2 = Array(cm.height()) { LongArray(cm.width()) }
            for (i in 0 until cm.height()) {
                for (j in 0 until cm.width())
                    gcm2[i][j] = gcm[p[i]][p[j]]
            }

            val awt = convertGCMtoAwt(IntegerMatrix(gcm2))
            val computer = IntegerMatrix.SNFComputer(awt)
            computer.doComputeSNF()
            val n = oddColumns.size
            val m = evenColumns.size
            print("\nExample #$name: ")
            println("Odd: $n; Even: $m;")

            println("\nSmith normal form: \n" + Arrays.toString(computer.s.getDiagonalEntries()))

            val nn = oddColumns.size + evenColumns.size

            assert(nn == Math.min(awt.width(), awt.height()))

            for (i in 0 until nn) {
                val entry = Math.abs(computer.s.myCoo[i][i])
                val parity = computer.isEven[i]
                if (parity) {
                    if (entry == 0L) {
                        print("K_2^{MW}(F); ")
                    } else if (entry == 2L) {
                        print("K_2^{MW}/<{u^2, v}>; ")
                    } else {
                        print("K_2^{MW}/" + entry / 2 + "<{u^2, v}>; ")
                    }
                } else {
                    if (entry == 0L) {
                        print("K_2(F);")
                    } else if (entry != 1L) {
                        print("K_2(F) / " + entry + "K_2(F); ")
                    }
                }
            }
            println()

            println("\nQ matrix: \n" + computer.q)
            println("\nPA matrix: \n" + computer.p.mul(awt))
            println("Indices: " + Arrays.toString(p))
        }

        fun doAllGCMExperiments() {
            gcmexperiment("27", arrayOf(longArrayOf(2, -1, 0),
                                        longArrayOf(-3, 2, -2),
                                        longArrayOf(0, -1, 2)))

            gcmexperiment("40", arrayOf(longArrayOf(2, -1, -2),
                                        longArrayOf(-1, 2, -2),
                                        longArrayOf(-2, -2, 2)))

            gcmexperiment("55", arrayOf(longArrayOf(2, -1, -2),
                                        longArrayOf(-2, 2, -2),
                                        longArrayOf(-2, -1, 2)))

            gcmexperiment("59", arrayOf(longArrayOf(2, -1, -2),
                                        longArrayOf(-2, 2, -2),
                                        longArrayOf(-2, -2, 2)))
            gcmexperiment("63", arrayOf(longArrayOf(2, -2, -3),
                                        longArrayOf(-1, 2, -2),
                                        longArrayOf(-1, -2, 2)))

            gcmexperiment("67", arrayOf(longArrayOf(2, -2, -4),
                                        longArrayOf(-1, 2, -2),
                                        longArrayOf(-1, -2, 2)))

            gcmexperiment("81", arrayOf(longArrayOf(2, -2, -3),
                                        longArrayOf(-2, 2, -2),
                                        longArrayOf(-1, -2, 2)))

            gcmexperiment("82", arrayOf(longArrayOf(2, -2, -4),
                                        longArrayOf(-2, 2, -2),
                                        longArrayOf(-1, -2, 2)))

            gcmexperiment("88", arrayOf(longArrayOf(2, -2, -1),
                                        longArrayOf(-2, 2, -1),
                                        longArrayOf(-4, -3, 2)))

            gcmexperiment("90", arrayOf(longArrayOf(2, -1, -2),
                                        longArrayOf(-4, 2, -4),
                                        longArrayOf(-2, -1, 2)))

            gcmexperiment("103", arrayOf(longArrayOf(2, -2, 0),
                                         longArrayOf(-2, 2, -1),
                                         longArrayOf(0, -1, 2)))

            gcmexperiment("106", arrayOf(longArrayOf(2, -2, 0),
                                         longArrayOf(-2, 2, -2),
                                         longArrayOf(0, -1, 2)))

            gcmexperiment("110", arrayOf(longArrayOf(2, -1, 0),
                                         longArrayOf(-4, 2, -2),
                                         longArrayOf(0, -1, 2)))

            gcmexperiment("113", arrayOf(longArrayOf(2, -2, 0),
                                         longArrayOf(-2, 2, -1),
                                         longArrayOf(0, -3, 2)))

            gcmexperiment("114", arrayOf(longArrayOf(2, -2, 0),
                                         longArrayOf(-2, 2, -3),
                                         longArrayOf(0, -1, 2)))

            gcmexperiment("116", arrayOf(longArrayOf(2, -1, 0),
                                         longArrayOf(-4, 2, -2),
                                         longArrayOf(0, -2, 2)))

            gcmexperiment("119", arrayOf(longArrayOf(2, -1, 0),
                                         longArrayOf(-3, 2, -4),
                                         longArrayOf(0, -1, 2)))

            gcmexperiment("122", arrayOf(longArrayOf(2, -1, 0),
                                         longArrayOf(-4, 2, -4),
                                         longArrayOf(0, -1, 2)))


            gcmexperiment("128", arrayOf(longArrayOf(2, -1, -1, 0),
                                         longArrayOf(-1, 2, -1, 0),
                                         longArrayOf(-1, -1, 2, -2),
                                         longArrayOf(0, 0, -1, 2)))

            gcmexperiment("135", arrayOf(longArrayOf(2, -1, 0, -1),
                                         longArrayOf(-1, 2, -2, 0),
                                         longArrayOf(0, -1, 2, -1),
                                         longArrayOf(-1, 0, -2, 2)))

            gcmexperiment("143", arrayOf(longArrayOf(2, -1, 0, -2),
                                         longArrayOf(-1, 2, -1, 0),
                                         longArrayOf(0, -2, 2, -2),
                                         longArrayOf(-1, 0, -1, 2)))

            gcmexperiment("144", arrayOf(longArrayOf(2, -1, 0, -2),
                                         longArrayOf(-1, 2, -2, 0),
                                         longArrayOf(0, -1, 2, -2),
                                         longArrayOf(-1, 0, -1, 2)))

            gcmexperiment("146", arrayOf(longArrayOf(2, -1, 0, -2),
                                         longArrayOf(-2, 2, -2, 0),
                                         longArrayOf(0, -1, 2, -2),
                                         longArrayOf(-1, 0, -1, 2)))

            gcmexperiment("147", arrayOf(longArrayOf(2, -1, 0, -2),
                                         longArrayOf(-2, 2, -1, 0),
                                         longArrayOf(0, -2, 2, -2),
                                         longArrayOf(-1, 0, -1, 2)))

            gcmexperiment("148", arrayOf(longArrayOf(2, -2, 0, -2),
                                         longArrayOf(-1, 2, -1, 0),
                                         longArrayOf(0, -2, 2, -2),
                                         longArrayOf(-1, 0, -1, 2)))

            gcmexperiment("156", arrayOf(longArrayOf(2, 0, 0, -1),
                                         longArrayOf(0, 2, 0, -1),
                                         longArrayOf(0, 0, 2, -1),
                                         longArrayOf(-2, -2, -2, 2)))

            gcmexperiment("162", arrayOf(longArrayOf(2, -1, 0, 0),
                                         longArrayOf(-1, 2, -1, 0),
                                         longArrayOf(0, -2, 2, -2),
                                         longArrayOf(0, 0, -1, 2)))

            gcmexperiment("168", arrayOf(longArrayOf(2, -1, 0, 0),
                                         longArrayOf(-2, 2, -1, 0),
                                         longArrayOf(0, -1, 2, -3),
                                         longArrayOf(0, 0, -1, 2)))

            gcmexperiment("174", arrayOf(longArrayOf(2, -1, 0, 0),
                                         longArrayOf(-2, 2, -2, 0),
                                         longArrayOf(0, -1, 2, -2),
                                         longArrayOf(0, 0, -1, 2)))

            gcmexperiment("179", arrayOf(longArrayOf(2, -1, -1, 0, 0),
                                         longArrayOf(-1, 2, 0, -1, 0),
                                         longArrayOf(-1, 0, 2, -1, 0),
                                         longArrayOf(0, -1, -1, 2, -2),
                                         longArrayOf(0, 0, 0, -1, 2)))

            gcmexperiment("197", arrayOf(longArrayOf(2, -1, 0, 0, 0),
                                         longArrayOf(-1, 2, -1, 0, 0),
                                         longArrayOf(0, -2, 2, -1, 0),
                                         longArrayOf(0, 0, -1, 2, -2),
                                         longArrayOf(0, 0, 0, -1, 2)))

            gcmexperiment("215", arrayOf(longArrayOf(2, -1, 0, 0, 0, 0),
                                         longArrayOf(-1, 2, -1, 0, 0, 0),
                                         longArrayOf(0, -2, 2, -1, 0, 0),
                                         longArrayOf(0, 0, -1, 2, -1, 0),
                                         longArrayOf(0, 0, 0, -1, 2, -2),
                                         longArrayOf(0, 0, 0, 0, -1, 2)))
        }
    }
}
