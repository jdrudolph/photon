using System;
using System.IO;
using System.Linq;
using System.Text;
using BaseLib.Graphic;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Generic;
using PerseusApi.Network;
using Path = System.IO.Path;

namespace PluginPHOTON
{
    public class PhotonFromNetwork : PluginInterop.Python.NetworkProcessing
    {
        public override string Heading => "Modifications";
        public override string Name => "PHOTON";
        public override string Description => "Reconstruct a signaling pathway";
        public override Bitmap2 DisplayImage => GraphUtils.ToBitmap2(Properties.Resources.icon);

        protected override string[] ReqiredPythonPackages => new[] { "perseuspy", "phos", "joblib" };
        public override int NumSupplTables => 2;
        public override DataType[] SupplDataTypes => new [] {DataType.Network, DataType.Matrix};

        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            codeFile = Path.GetTempFileName();
            File.WriteAllText(codeFile, Encoding.UTF8.GetString(Properties.Resources.photon_from_network));
            return true;
        }

        protected override Parameter[] SpecificParameters(INetworkData ndata, ref string errString)
        {
            if (ndata.Count() != 1)
            {
                errString = "Please make sure to have only one network in the collection.";
                return null;
            }
            var network = ndata.Single();
            var nodeTable = network.NodeTable;
            if (nodeTable.MultiNumericColumnCount < 1)
            {
                errString = "Please add a multi numeric column with the data to the node table";
                return null;
            }
            var edgeTable = network.EdgeTable;
            var confidenceColumn = edgeTable.NumericColumnNames.FindIndex(col => col.ToLower().Equals("confidence"));
            return new Parameter[]
            {
                new MultiChoiceParam("Data columns")
                {
                    Value = Enumerable.Range(0, nodeTable.MultiNumericColumnCount).Where(i => !new [] {"position"}.Contains(nodeTable.MultiNumericColumnNames[i].ToLower())).ToArray(),
                    Values = nodeTable.MultiNumericColumnNames,
                    Help = "PTM data of different conditions for the nodes in the network."
                },
                new SingleChoiceParam("Confidence column")
                {
                    Value = confidenceColumn == -1 ? edgeTable.NumericColumnCount : confidenceColumn,
                    Values = edgeTable.NumericColumnNames.Concat(new []{"Use constant value"}).ToList(),
                    Help = "Confidence score for interactions. Will be used as weights in the signaling score calculation"
                }, 
                new StringParam("Signaling source")
                {
                    Help = "Has to be a human entrez gene id. Select the starting point of the signaling network (optional)."
                }, 
                new IntParam("Required number of observations", 4) {Help = "Required minumum number of observations for score calculation."},
                new IntParam("Number of permutations", 1000) {Help = "Number of permutations used for empirical p-value calculation"},
            };
        }
    }
}