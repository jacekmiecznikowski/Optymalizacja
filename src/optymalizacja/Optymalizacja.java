package optymalizacja;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

class Population {

    public String binary_seq[][];
    public String long_binary[];
    public double double_seq[][];

    public Population(int k, int n) {
        binary_seq = new String[k][n];
        long_binary = new String[k];
        double_seq = new double[k][n];
    }

    public void print_binary_seq() {
        for (int i = 0; i < binary_seq.length; i++) {
            for (int j = 0; j < binary_seq[0].length; j++) {
                if (j == 0) {
                    System.out.print("x" + i + " = (" + binary_seq[i][j] + ", ");
                } else if (j == binary_seq[0].length - 1) {
                    System.out.print(binary_seq[i][j] + ")\n");
                } else {
                    System.out.print(binary_seq[i][j] + ", ");
                }
            }
        }
    }

    public void print_long_binary() {
        for (int i = 0; i < long_binary.length; i++) {
            System.out.println("x" + i + " = " + long_binary[i]);
        }
    }

    public void print_double() {
        for (int i = 0; i < double_seq.length; i++) {
            for (int j = 0; j < double_seq[0].length; j++) {
                if (j == 0) {
                    System.out.print("x" + i + " = (" + double_seq[i][j] + ", ");
                } else if (j == double_seq[0].length - 1) {
                    System.out.print(double_seq[i][j] + ")\n");
                } else {
                    System.out.print(double_seq[i][j] + ", ");
                }
            }
        }
    }
}

public class Optymalizacja {

    static double fun(double[][] x, int indeks, int n, double A, double w) {
        double val = A * n;
        for (int i = 0; i < n; i++) {
            val += Math.pow(x[indeks][i], 2) - A * Math.cos(w * x[indeks][i]);
        }
        //return val;
        return -val+100;
    }

    static void print_fun(Population population, int k, int n, double A, double w) {
        double sum = 0;
        for (int i = 0; i < k; i++) {
            System.out.println("f(x" + i + ") = " + fun(population.double_seq, i, n, A, w));
            sum += fun(population.double_seq, i, n, A, w);
        }
        System.out.println("F srednie " + sum/k);
    }
    static double avg(Population population, int k, int n, double A, double w) {
        double sum = 0;
        for (int i = 0; i < k; i++) {
            sum += fun(population.double_seq, i, n, A, w);
        }
       return sum/k;
    }

    static String generate_binary(int len) {
        String binary = "";
        for (int i = 0; i < len; i++) {
            binary += 0 + (int) (Math.random() * ((1 - 0) + 1));
        }
        return binary;
    }

    static int binary_to_int(String binary) {
        char[] numbers = binary.toCharArray();
        int integer = 0;
        for (int i = numbers.length - 1; i >= 0; i--) {
            if (numbers[i] == '1') {
                integer += Math.pow(2, (numbers.length - i - 1));
            }
        }
        return integer;
    }

    static double translate(int num, double a, double b, int m) {
        return a + (b - a) * num / (Math.pow(2, m) - 1);
    }

    static String[][] longBinaryToBinarySeq(String[] long_binary, int[] m, int n) {
        String[][] binary_seq = new String[long_binary.length][n];
        for (int i = 0; i < long_binary.length; i++) {
            String tmp_binary = long_binary[i];
            for (int j = 0; j < n; j++) {
                binary_seq[i][j] = tmp_binary.substring(0, m[j]);
                tmp_binary = tmp_binary.substring(m[j]);
            }

        }
        return binary_seq;
    }

    static String[] binarySeqToLongBinary(String[][] binary_seq) {
        String[] long_binary = new String[binary_seq.length];
        for (int i = 0; i < binary_seq.length; i++) {
            String tmp = "";
            for (String binary_seq1 : binary_seq[i]) {
                tmp += binary_seq1;
            }
            long_binary[i] = tmp;

        }
        return long_binary;

    }

    static boolean contains(int[] array, int v) {
        for (int e : array) {
            if (e == v) {
                return true;
            }
        }

        return false;
    }

    static String replace(String binary_seq, int index, char replace) {
        if (binary_seq == null) {
            return binary_seq;
        } else if (index < 0 || index >= binary_seq.length()) {
            return binary_seq;
        }
        char[] chars = binary_seq.toCharArray();
        chars[index] = replace;
        return String.valueOf(chars);
    }

    static Population ranking(Population population, int n, double A, double w) {
        int smaller = population.double_seq.length - 1;
        int bigger = smaller - 1;
        int index;
        double tmp_double;
        String tmp_binary;
        String tmp_long_binary;
        while (bigger >= 0) {
            if (fun(population.double_seq, bigger, n, A, w) < fun(population.double_seq, smaller, n, A, w)) {
                for (int i = 0; i < n; i++) {
                    tmp_double = population.double_seq[bigger][i];
                    tmp_binary = population.binary_seq[bigger][i];
                    population.double_seq[bigger][i] = population.double_seq[smaller][i];
                    population.binary_seq[bigger][i] = population.binary_seq[smaller][i];
                    population.double_seq[smaller][i] = tmp_double;
                    population.binary_seq[smaller][i] = tmp_binary;
                }
                tmp_long_binary = population.long_binary[bigger];
                population.long_binary[bigger] = population.long_binary[smaller];
                population.long_binary[smaller] = tmp_long_binary;
                smaller = population.double_seq.length - 1;
                bigger = smaller - 1;
            } else {
                smaller--;
                bigger--;
            }
        }
        Population new_population = new Population(population.double_seq.length, n);
        for (int i = 0; i < population.double_seq.length; i++) {
            index = (0 + (int) (Math.random() * (((population.double_seq.length - 1) - 0) + 1)));
            index = (0 + (int) (Math.random() * (((index) - 0) + 1)));
            for (int j = 0; j < population.double_seq[0].length; j++) {
                new_population.double_seq[i][j] = population.double_seq[index][j];
                new_population.binary_seq[i][j] = population.binary_seq[index][j];
            }
            new_population.long_binary[i] = population.long_binary[index];
        }
        return new_population;
    }

    static Population tournament(Population population, int reuse, int n, double A, double w) {
        int group_size = (int) (0.2 * population.binary_seq.length);
        Population new_population = new Population(population.double_seq.length, n);
        Population group = new Population(group_size, n);
        for (int i = 0; i < population.double_seq.length; i++) {
            int[] group_indexes = new int[group_size];
            for (int j = 0; j < group.double_seq.length; j++) {
                if (reuse == 1) {
                    int random = (0 + (int) (Math.random() * (((population.double_seq.length - 1) - 0) + 1)));
                    for (int k = 0; k < group.double_seq[0].length; k++) {
                        group.double_seq[j][k] = population.double_seq[random][k];
                        group.binary_seq[j][k] = population.binary_seq[random][k];
                    }
                    group.long_binary[j] = population.long_binary[random];

                } else {
                    int random = (0 + (int) (Math.random() * (((population.double_seq.length - 1) - 0) + 1)));
                    for (int k = 0; k < group.double_seq[0].length; k++) {
                        do {
                            random = (0 + (int) (Math.random() * (((population.double_seq.length - 1) - 0) + 1)));
                        } while (contains(group_indexes, random));
                        group_indexes[k] = random;
                        group.double_seq[j][k] = population.double_seq[random][k];
                        group.binary_seq[j][k] = population.binary_seq[random][k];
                    }
                    group.long_binary[j] = population.long_binary[random];
                }
            }
            int max = 0;
            for (int j = 1; j < group.double_seq.length; j++) {
                double newnumber = fun(group.double_seq, j, n, A, w);
                if (newnumber > fun(group.double_seq, max, n, A, w)) {
                    max = j;
                }
            }
            for (int j = 0; j < population.double_seq[0].length; j++) {
                new_population.double_seq[i][j] = group.double_seq[max][j];
                new_population.binary_seq[i][j] = group.binary_seq[max][j];
            }
            new_population.long_binary[i] = group.long_binary[max];
        }
        return new_population;
    }

    static Population roulette(Population population, int n, double A, double w) {
        double F = 0;
        double p[] = new double[population.double_seq.length];
        double q[] = new double[population.double_seq.length];
        for (int i = 0; i < population.double_seq.length; i++) {
            F += fun(population.double_seq, i, n, A, w);
        }
        for (int i = 0; i < population.double_seq.length; i++) {
            p[i] = fun(population.double_seq, i, n, A, w) / F;
        }
        for (int i = 0; i < population.double_seq.length; i++) {
            q[i] = 0;
            for (int j = 0; j <= i; j++) {
                q[i] += p[j];
            }
        }
        Population new_population = new Population(population.double_seq.length, n);
        for (int i = 0; i < population.double_seq.length; i++) {
            double r = 0 + (1 - 0) * Math.random();
            int index = 0;
            do {
                if (r < q[index]) {
                    for (int j = 0; j < population.double_seq[0].length; j++) {
                        new_population.double_seq[i][j] = population.double_seq[index][j];
                        new_population.binary_seq[i][j] = population.binary_seq[index][j];
                    }
                    new_population.long_binary[i] = population.long_binary[index];
                    break;
                }
                index++;
            } while (true);

        }

        return new_population;
    }

    static Population mutation(double pm, Population population, int n, double A, double w, double[] a, double[] b, int[] m) {
        Population new_population = new Population(population.double_seq.length, n);
        double r;
        for (int i = 0; i < population.binary_seq.length; i++) {
            for (int j = 0; j < population.binary_seq[0].length; j++) {
                new_population.binary_seq[i][j] = population.binary_seq[i][j];
                for (int k = 0; k < population.binary_seq[i][j].length(); k++) {
                    r = Math.random();
                    if (r < pm) {
                        new_population.binary_seq[i][j] = replace(new_population.binary_seq[i][j], k, (population.binary_seq[i][j].charAt(k) == 0) ? '1' : '0');
                    }
                }
                new_population.double_seq[i][j] = translate(binary_to_int(new_population.binary_seq[i][j]), a[j], b[j], m[j]);
            }
        }
        new_population.long_binary = binarySeqToLongBinary(new_population.binary_seq);
        return new_population;
    }

    static Population inversion(double pi, Population population, int n, double A, double w, double[] a, double[] b, int[] m) {
        Population new_population = new Population(population.double_seq.length, n);
        for (int i = 0; i < population.binary_seq.length; i++) {
            for (int j = 0; j < population.binary_seq[0].length; j++) {
                new_population.binary_seq[i][j] = population.binary_seq[i][j];
                double r = 0 + (1 - 0) * Math.random();
                if (r < pi) {
                    int from = 0 + (int) (Math.random() * ((m[j] - 1) + 1));
                    int to, temp;
                    do {
                        to = 0 + (int) (Math.random() * ((m[j] - 1) + 1));
                    } while (to == from);
                    if (to < from) {
                        temp = from;
                        from = to;
                        to = temp;
                    }
                    StringBuilder dest = new StringBuilder(population.binary_seq[i][j].length());
                    for (int k = to; k >= from; k--) {
                        dest.append(population.binary_seq[i][j].charAt(k));
                    }
                    if (from > 0) {
                        dest.insert(0, population.binary_seq[i][j].substring(0, from));
                        dest.append(population.binary_seq[i][j].substring(to + 1));
                    } else {
                        dest.append(population.binary_seq[i][j].substring(to + 1));
                    }
                    new_population.binary_seq[i][j] = dest.toString();
                }
                new_population.double_seq[i][j] = translate(binary_to_int(new_population.binary_seq[i][j]), a[j], b[j], m[j]);
            }
        }
        new_population.long_binary = binarySeqToLongBinary(new_population.binary_seq);
        return new_population;
    }

    static Population crossover(double pc, int pointsNum, Population population, int n, double A, double w, int sum, double[] a, double[] b, int[] m) {
        Population new_population = new Population(population.double_seq.length, n);
        List indexList = new ArrayList();
        for (int i = 0; i < population.long_binary.length; i++) {
            double r = 0 + (1 - 0) * Math.random();
            if (r < pc) {
                indexList.add(i);
            }
        }
        Collections.shuffle(indexList);
        if (indexList.size() % 2 != 0) {
            indexList.remove(0 + (int) (Math.random() * ((indexList.size() - 1) + 1)));
        }
        System.out.println("Krzyżowane pary: " + indexList);
        if (pointsNum > 0) {
            int[] points = new int[pointsNum];
            int r;
            System.out.println("Punkt(y) przecinania:");
            for (int i = 0; i < pointsNum; i++) {
                do {
                    r = 0 + (int) (Math.random() * (sum));
                } while (contains(points, r));
                points[i] = r;
                System.out.println(points[i]);
            }
            Arrays.sort(points);

            for (int i = 0; i < indexList.size(); i += 2) {
                new_population.long_binary[(int) indexList.get(i)] = "";
                new_population.long_binary[(int) indexList.get(i + 1)] = "";
                for (int j = 0; j <= pointsNum; j++) {
                    if (j == 0) {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i)].substring(0, points[0] + 1);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i + 1)].substring(0, points[0] + 1);
                    } else if ((j == pointsNum) && (j % 2 == 0)) {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i)].substring(points[j - 1] + 1);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i + 1)].substring(points[j - 1] + 1);
                    } else if ((j == pointsNum) && (j % 2 != 0)) {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i + 1)].substring(points[j - 1] + 1);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i)].substring(points[j - 1] + 1);
                    } else if (j % 2 == 0) {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i)].substring(points[j - 1] + 1, points[j] + 1);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i + 1)].substring(points[j - 1] + 1, points[j] + 1);
                    } else if (j % 2 != 0) {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i + 1)].substring(points[j - 1] + 1, points[j] + 1);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i)].substring(points[j - 1] + 1, points[j] + 1);
                    }
                }

                population.long_binary[(int) indexList.get(i)] = new_population.long_binary[(int) indexList.get(i)];
                population.long_binary[(int) indexList.get(i + 1)] = new_population.long_binary[(int) indexList.get(i + 1)];

            }

        } else {
            String pattern = generate_binary(sum);
            System.out.println("Wzór krzyżowania: " + pattern);
            char[] patternArray = pattern.toCharArray();

            for (int i = 0; i < indexList.size(); i += 2) {
                new_population.long_binary[(int) indexList.get(i)] = "";
                new_population.long_binary[(int) indexList.get(i + 1)] = "";
                for (int j = 0; j < patternArray.length; j++) {
                    if (patternArray[j] == '0') {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i)].charAt(j);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i + 1)].charAt(j);
                    } else {
                        new_population.long_binary[(int) indexList.get(i)] += population.long_binary[(int) indexList.get(i + 1)].charAt(j);
                        new_population.long_binary[(int) indexList.get(i + 1)] += population.long_binary[(int) indexList.get(i)].charAt(j);
                    }
                }
                population.long_binary[(int) indexList.get(i)] = new_population.long_binary[(int) indexList.get(i)];
                population.long_binary[(int) indexList.get(i + 1)] = new_population.long_binary[(int) indexList.get(i + 1)];
            }
        }
        population.binary_seq = longBinaryToBinarySeq(population.long_binary, m, n);
        for (int i = 0; i < population.binary_seq.length; i++) {
            for (int j = 0; j < population.binary_seq[i].length; j++) {
                population.double_seq[i][j] = translate(binary_to_int(population.binary_seq[i][j]), a[j], b[j], m[j]);
            }
        }
        return population;
    }

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        double A, w;
        int n = 2;
        int k = 15;
        int generations = 200;
        double[] a = {-4.0, -4.0};
        double[] b = {4.0, 4.0};
        double[] d = {6.0, 6.0};
        int[] m = new int[n];
        A = 10;
        w = 2 * Math.PI;
        int sum = 0;
        double value;
        double pi=0.1;
        double pm = 0.1;
        double pc= 0.5;
        int cross= 0;

        for (int i = 0; i < n; i++) {
            value = (b[i] - a[i]) * Math.pow(10, d[i]);
            m[i] = 0;
            while (value >= Math.pow(2, m[i])) {
                m[i]++;
            }
            sum += m[i];
        }
        System.out.println("Łączna długość łańcuchów m[] =  " + sum);
        double avg = 0;
        for (int test=0; test<5;test++) {
        Population population = new Population(k, n);
        Population new_population = new Population(k, n);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                population.binary_seq[i][j] = generate_binary(m[j]);
                population.double_seq[i][j] = translate(binary_to_int(population.binary_seq[i][j]), a[j], b[j], m[j]);

            }
        }
        population.long_binary = binarySeqToLongBinary(population.binary_seq);
        System.out.println("Wygenerowana populacja:");
        population.print_binary_seq();
        population.print_long_binary();
        population.print_double();
        print_fun(population, k, n, A, w);
        System.out.println("Wybierz metodę selekcji.");
        System.out.println("1. Ranking\n" + "2. Turniej\n" + "3. Ruletka");
       String method ="1";// String method = in.next();
        String choice = "";
        System.out.println("Wybierz metodę sukcesji.");
        System.out.println("1. Z całkowitym zastępowaniem\n" + "2. Losowa\n" + "3. Elitarna");
       String type= "3"; //String type = in.next();
        switch (type) {
            case "1":
                for (int i = 0; i < generations; i++) {

                    switch (method) {
                        case "1":
                            population = ranking(population, n, A, w);
                            break;
                        case "2":
                            if (choice.equals("")) {
                                System.out.println("Ze zwracaniem?\n" + "1. Tak\n" + "2. Nie");
                                choice = in.next();
                            }
                            switch (choice) {
                                case "1":
                                case "2":
                                    population = tournament(population, Integer.parseInt(method), n, A, w);
                                    break;
                                default:
                                    System.out.println("Błąd.");
                                    System.exit(0);
                            }
                            break;
                        case "3":
                            population = roulette(population, n, A, w);
                            break;
                        default:
                            System.out.println("Błąd.");
                            System.exit(0);
                    }
                    System.out.println("Po selekcji:");
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);

                    population = mutation(pm, population, n, A, w, a, b, m);
                    System.out.println("Po mutacji:");
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);

                    population = inversion(pi, population, n, A, w, a, b, m);
                    System.out.println("Po inwersji:");
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);

                    System.out.println("Po krzyżowaniu:");
                    population = crossover(pc, cross, population, n, A, w, sum, a, b, m);
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);

                }
                break;
            case "2":
                for (int i = 0; i < generations; i++) {

                    switch (method) {
                        case "1":
                            new_population = ranking(population, n, A, w);
                            break;
                        case "2":
                            if (choice.equals("")) {
                                System.out.println("Ze zwracaniem?\n" + "1. Tak\n" + "2. Nie");
                                choice = in.next();
                            }
                            switch (choice) {
                                case "1":
                                case "2":
                                    new_population = tournament(population, Integer.parseInt(method), n, A, w);
                                    break;
                                default:
                                    System.out.println("Błąd.");
                                    System.exit(0);
                            }
                            break;
                        case "3":
                            new_population = roulette(population, n, A, w);
                            break;
                        default:
                            System.out.println("Błąd.");
                            System.exit(0);
                    }
                    System.out.println("Po selekcji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    new_population = mutation(pm, new_population, n, A, w, a, b, m);
                    System.out.println("Po mutacji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    new_population = inversion(pi, new_population, n, A, w, a, b, m);
                    System.out.println("Po inwersji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    System.out.println("Po krzyżowaniu:");
                    new_population = crossover(pc, cross, new_population, n, A, w, sum, a, b, m);
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);
                    for (int j = 0; j < k; j++) {
                        if (0 + (int) (Math.random() * ((1 - 0) + 1)) == 0) {
                            population.long_binary[j] = new_population.long_binary[j];
                            for (int l = 0; l < population.binary_seq[j].length; l++) {
                                population.binary_seq[j][l] = new_population.binary_seq[j][l];
                                population.double_seq[j][l] = new_population.double_seq[j][l];
                            }
                        }
                    }
                    System.out.println("Po sukcesji:");
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);
                }
                break;
            case "3":
                for (int i = 0; i < generations; i++) {

                    switch (method) {
                        case "1":
                            new_population = ranking(population, n, A, w);
                            break;
                        case "2":
                            if (choice.equals("")) {
                                System.out.println("Ze zwracaniem?\n" + "1. Tak\n" + "2. Nie");
                                choice = in.next();
                            }
                            switch (choice) {
                                case "1":
                                case "2":
                                    new_population = tournament(population, Integer.parseInt(method), n, A, w);
                                    break;
                                default:
                                    System.out.println("Błąd.");
                                    System.exit(0);
                            }
                            break;
                        case "3":
                            new_population = roulette(population, n, A, w);
                            break;
                        default:
                            System.out.println("Błąd.");
                            System.exit(0);
                    }
                    System.out.println("Po selekcji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    new_population = mutation(pm, new_population, n, A, w, a, b, m);
                    System.out.println("Po mutacji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    new_population = inversion(pi, new_population, n, A, w, a, b, m);
                    System.out.println("Po inwersji:");
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    System.out.println("Po krzyżowaniu:");
                    new_population = crossover(pc, cross, new_population, n, A, w, sum, a, b, m);
                    new_population.print_binary_seq();
                    new_population.print_long_binary();
                    new_population.print_double();
                    print_fun(new_population, k, n, A, w);

                    Population merged_population = new Population(k * 2, n);

                    for (int j = 0; j < k * 2; j++) {
                        for (int l = 0; l < population.binary_seq[0].length; l++) {
                            if (j < k) {
                                merged_population.binary_seq[j][l] = new_population.binary_seq[j][l];
                                merged_population.double_seq[j][l] = new_population.double_seq[j][l];
                            } else {
                                merged_population.binary_seq[j][l] = new_population.binary_seq[j - k][l];
                                merged_population.double_seq[j][l] = new_population.double_seq[j - k][l];
                            }
                        }
                        if (j < k) {
                            merged_population.long_binary[j] = new_population.long_binary[j];
                        } else {
                            merged_population.long_binary[j] = new_population.long_binary[j - k];
                        }
                    }
                    int smaller = merged_population.double_seq.length - 1;
                    int bigger = smaller - 1;
                    double tmp_double;
                    String tmp_binary;
                    String tmp_long_binary;
                    while (bigger >= 0) {
                        if (fun(merged_population.double_seq, bigger, n, A, w) < fun(merged_population.double_seq, smaller, n, A, w)) {
                            for (int z = 0; z < n; z++) {
                                tmp_double = merged_population.double_seq[bigger][z];
                                tmp_binary = merged_population.binary_seq[bigger][z];
                                merged_population.double_seq[bigger][z] = merged_population.double_seq[smaller][z];
                                merged_population.binary_seq[bigger][z] = merged_population.binary_seq[smaller][z];
                                merged_population.double_seq[smaller][z] = tmp_double;
                                merged_population.binary_seq[smaller][z] = tmp_binary;
                            }
                            tmp_long_binary = merged_population.long_binary[bigger];
                            merged_population.long_binary[bigger] = merged_population.long_binary[smaller];
                            merged_population.long_binary[smaller] = tmp_long_binary;
                            smaller = merged_population.double_seq.length - 1;
                            bigger = smaller - 1;
                        } else {
                            smaller--;
                            bigger--;
                        }
                    }
                    for (int j = 0; j < k; j++) {
                        population.long_binary[j] = merged_population.long_binary[j];
                        for (int l = 0; l < population.binary_seq[j].length; l++) {
                            population.binary_seq[j][l] = merged_population.binary_seq[j][l];
                            population.double_seq[j][l] = merged_population.double_seq[j][l];
                        }
                    }
                    System.out.println("Po sukcesji:");
                    population.print_binary_seq();
                    population.print_long_binary();
                    population.print_double();
                    print_fun(population, k, n, A, w);
                }
                break;
            default:
                System.out.println("Błąd.");
                System.exit(0);
        }
        avg+=avg(population, k, n, A, w);
        }
        avg = avg/5;
        System.out.println("AVG: "+avg);
    }
}