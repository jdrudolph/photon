﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using BaseLib.Graphic;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Generic;
using PerseusApi.Network;
using PluginInterop;
using Path = System.IO.Path;

namespace PluginPHOTON
{
    public class AnatFromNetwork : PluginInterop.Python.NetworkProcessing
    {
        public override string Heading => "Subnetwork";
        public override string Name => "ANAT";
        public override string Description => "Subnetwork reconstruction using the Advanced Network Analysis Tool (ANAT).";
        public override Bitmap2 DisplayImage => null;

        protected override string[] ReqiredPythonPackages => new[] { "perseuspy", "phos" };
        public override int NumSupplTables => 1;
        public override DataType[] SupplDataTypes => new [] {DataType.Network};

        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            codeFile = Path.GetTempFileName();
            File.WriteAllText(codeFile, Encoding.UTF8.GetString(Properties.Resources.anat_from_network));
            return true;
        }

		protected override string GetCommandLineArguments(Parameters param)
		{
			var tempFile = Path.GetTempFileName();
			param.ToFile(tempFile);
			return tempFile;
		}

		protected override Parameter[] SpecificParameters(INetworkDataAnnColumns ndata, ref string errString)
        {
            if (ndata.Count() != 1)
            {
                errString = "Please make sure to have only one network in the collection.";
                return null;
            }
            var network = ndata.Single();
            var nodeTable = network.NodeTable;
            if (nodeTable.CategoryColumnCount < 1)
            {
                errString = "Please add at least one category column to the node table.";
                return null;
            }
            var edgeTable = network.EdgeTable;
            var confidenceColumn = edgeTable.NumericColumnNames.FindIndex(col => col.ToLower().Equals("confidence"));
	        var signalingSourceParam = new StringParam("Signaling source")
	        {
		        Help = "Select the starting point (anchor) of the signaling network. Leave blank for unrooted network. " +
		               "Has to be a Node from the network. Choose multiple anchors by separating ids by a space."
	        };
            var terminalsParam = new SingleChoiceWithSubParams("Signaling targets", 0)
            {
                Values = nodeTable.CategoryColumnNames,
                SubParams = Enumerable.Range(0, nodeTable.CategoryColumnCount).Select(i =>
                {
                    var param = new MultiChoiceParam(nodeTable.CategoryColumnNames[i])
                    {
                        Values = nodeTable.GetCategoryColumnValuesAt(i)
                    };
                    return new Parameters(param);
                }).ToArray(),
                Help = "Select the end points (terminals) of the signaling network. Using more than 100 terminals is not recommended.",
            };
            var confidenceColumnParam = new SingleChoiceParam("Confidence column")
            {
                Value = confidenceColumn == -1 ? edgeTable.NumericColumnCount : confidenceColumn,
                Values = edgeTable.NumericColumnNames.Concat(new []{"Use constant value"}).ToList(),
                Help = "Confidence score for interactions. Will be used as weights in the network reconstruction"
            };
            return new Parameter[]
            {
                confidenceColumnParam,
                signalingSourceParam,
                terminalsParam
            };
        }
    }
}