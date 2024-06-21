import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class FordFulkerson {
    public static void main(String[] args) {
        if (args.length != 1) {
            System.out.println("Usage: java minSPT <input_file>");
            return;
        }
        // Get the input file name from the command-line arguments
        String inputFileName = args[0];
        try {
            List<Graph> graphs = parseGraphs(inputFileName);

            for (int i = 0; i < graphs.size(); i++) {
                Graph graph = graphs.get(i);
                int source = 0; // Replace with the actual source vertex index
                int sink = graph.getSize() - 1; // Replace with the actual sink vertex index
                double maxFlow = fordFulkerson(graph, source, sink);

                // Print the graph details
                System.out.println("** G" + (i + 1) + ": |V|=" + graph.getSize());
                if (graph.getSize() <= 10) {
                    printFlowNetwork(graph);
                }

                // printFlowGraph(flowGraph, graph.getSize());
                long startTime = System.currentTimeMillis();
                // int maxFlow = algorithm.preflowPush(0, graph.getSize() - 1, graph.getSize());
                long endTime = System.currentTimeMillis();

                System.out.println("Max flow ==> " + maxFlow + "(" + (endTime - startTime) + " ms)");
                System.out.println(" ");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("*** Asg 7 by Drithi Madagani ***");
    }

    private static void printFlowNetwork(Graph graph) {
        System.out.println("       Flow network:");
        System.out.printf("%10s", " "); // Adjust the indentation
        for (int i = 0; i < graph.getSize(); i++) {
            System.out.printf("%4d", i);
        }
        System.out.println();
        System.out.print("         ");
        for (int i = 0; i < graph.getSize(); i++) {
            System.out.print("----");
        }
        System.out.println();
        for (int i = 0; i < graph.getSize(); i++) {
            System.out.printf("%4d:  ", i);
            List<Edge> edges = graph.getEdges();
            for (int j = 0; j < graph.getSize(); j++) {
                double weight = 0;
                for (Edge edge : edges) {
                    if (edge.u == i && edge.v == j) {
                        weight = edge.weight;
                        break;
                    }
                }
                System.out.printf("%4s", weight != 0 ? String.format("%4.0f", weight) : "-");
            }
            System.out.println();
        }
    }

    private static void printFlowGraph(double[][] flowGraph, int size) {
        System.out.println("Maximum Flow :");
        for (int u = 0; u < size; u++) {
            for (int v = 0; v < size; v++) {
                System.out.print(flowGraph[u][v] + " ");
            }
            System.out.println();
        }
    }

    public static double fordFulkerson(Graph graph, int source, int sink) {
        // Basic input validation
        if (graph == null || source < 0 || sink < 0 || source >= graph.getSize() || sink >= graph.getSize()) {
            throw new IllegalArgumentException("Invalid graph or source/sink nodes");
        }

        Map<Integer, Map<Integer, Double>> residualMap = new HashMap<>(); // Map for residual graph
        Map<Integer, Map<Integer, Double>> flowMap = new HashMap<>(); // Map for flow graph
        double maxFlow = 0; // To store the final flow value

        // Initialize the residual graph with the capacities from the input graph
        graph.getEdges().forEach(edge -> {
            residualMap.computeIfAbsent(edge.u, k -> new HashMap<>()).put(edge.v, edge.weight);
            flowMap.computeIfAbsent(edge.u, k -> new HashMap<>()).putIfAbsent(edge.v, 0.0);
        });

        int[] parent = new int[graph.getSize()]; // Array to store the augmenting path

        // Augment the flow as long as there is a path from source to sink
        while (bfs(residualMap, source, sink, parent)) {
            double pathFlow = Double.POSITIVE_INFINITY; // Find the maximum flow through the path found.

            // Find the minimum residual capacity of the edges along the path
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = Math.min(pathFlow, residualMap.get(u).get(v));
            }

            // Update the residual capacities and add to the flow graph
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                residualMap.get(u).put(v, residualMap.get(u).get(v) - pathFlow);
                residualMap.computeIfAbsent(v, k -> new HashMap<>()).merge(u, pathFlow, Double::sum);
                flowMap.get(u).merge(v, pathFlow, Double::sum);
            }

            maxFlow += pathFlow; // Add the path flow to the total flow
        }

        return maxFlow; // Return the total flow
    }

    private static void printFlowGraph(double[][] flowGraph) {
        for (int i = 0; i < flowGraph.length; i++) {
            for (int j = 0; j < flowGraph[i].length; j++) {
                if (flowGraph[i][j] > 0) {
                    System.out.println("Flow from " + i + " to " + j + ": " + flowGraph[i][j]);
                }
            }
        }
    }

    private static boolean bfs(Map<Integer, Map<Integer, Double>> residualGraph, int source, int sink, int[] parent) {
        BitSet visited = new BitSet();
        Deque<Integer> deque = new ArrayDeque<>();

        Arrays.fill(parent, -1);
        deque.offer(source);
        visited.set(source);

        while (!deque.isEmpty()) {
            int u = deque.poll();
            Map<Integer, Double> edges = residualGraph.get(u);

            if (edges != null) {
                for (Map.Entry<Integer, Double> entry : edges.entrySet()) {
                    int v = entry.getKey();
                    double capacity = entry.getValue();

                    if (!visited.get(v) && capacity > 0) {
                        deque.offer(v);
                        parent[v] = u;
                        visited.set(v);

                        if (v == sink) { // Early exit if sink is found
                            return true;
                        }
                    }
                }
            }
        }

        return visited.get(sink);
    }

    private static List<Graph> parseGraphs(String path) throws IOException {
        List<Graph> graphs = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            boolean inGraphSection = false;
            Graph currentGraph = null;

            while ((line = reader.readLine()) != null) {
                if (line.matches("\\*\\* G\\d+:.*")) {
                    if (currentGraph != null) {
                        graphs.add(currentGraph);
                    }
                    currentGraph = new Graph();
                    inGraphSection = true;
                } else if (inGraphSection && !line.startsWith("-")) {
                    parseGraphEdges(currentGraph, line);
                } else if (inGraphSection && line.startsWith("-")) {
                    inGraphSection = false;
                }
            }

            if (currentGraph != null) {
                graphs.add(currentGraph);
            }
        }
        return graphs;
    }

    // 11/12/2023 to 11/15/2023
    // Function to parse edges and add them to a graph
    private static void parseGraphEdges(Graph graph, String line) {
        if (line.contains("(")) {
            String[] parts = line.substring(line.indexOf('(') + 1, line.indexOf(')')).split(",");
            try {
                int u = Integer.parseInt(parts[0].trim());
                int v = Integer.parseInt(parts[1].trim());
                double weight = Double.parseDouble(parts[2].trim());
                graph.addEdge(u, v, weight);
            } catch (NumberFormatException e) {
                // Handle parsing errors
            }
        }
    }
}

class Graph {
    private List<Edge> edges;
    private int size;

    // Constructor to initialize the graph
    public Graph() {
        size = 0;
        edges = new ArrayList<>();
    }

    // Method to add an edge to the graph
    public void addEdge(int u, int v, double weight) {
        edges.add(new Edge(u, v, weight));
        size = Math.max(size, Math.max(u, v) + 1);
    }

    public int getSize() {
        return size;
    }

    public List<Edge> getEdges() {
        return edges;
    }
}

class Edge implements Comparable<Edge> {
    int u;
    int v;
    double weight;
    double flow;

    Edge(int u, int v, double weight) {
        this.u = u;
        this.v = v;
        this.weight = weight;
        this.flow = 0;
    }

    @Override
    public int compareTo(Edge other) {
        return Double.compare(this.weight, other.weight);
    }
}
