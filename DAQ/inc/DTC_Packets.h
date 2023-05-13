#ifndef DTC_PACKETS_H
#define DTC_PACKETS_H

#include <bitset>
#include <cstdint>  // uint8_t, uint16_t
#include <vector>
#include <cassert>

#include "mu2e_pcie_utils/dtcInterfaceLib/DTC_Types.h"

#include "mu2e_driver/mu2e_mmap_ioctl.h"

namespace DTCLib {

/// <summary>
/// Defined Packet Types for the DTC DMA Protocol
/// </summary>
enum DTC_PacketType : uint8_t
{
	DTC_PacketType_DCSRequest = 0,
	DTC_PacketType_Heartbeat = 1,
	DTC_PacketType_DataRequest = 2,
	DTC_PacketType_DCSReply = 4,
	DTC_PacketType_DataHeader = 5,
	DTC_PacketType_Invalid = 0x10,
};

/// <summary>
/// Possible values for the Status word of the Data Header packet
/// </summary>
enum DTC_DataStatus
{
	DTC_DataStatus_Valid = 0,
	DTC_DataStatus_NoValid = 1,
	DTC_DataStatus_Invalid = 2,
};

/// <summary>
/// Possible values for the Op word of the DCS Request packet.
/// </summary>
enum DTC_DCSOperationType : uint8_t
{
	DTC_DCSOperationType_Read = 0,
	DTC_DCSOperationType_Write = 1,
	DTC_DCSOperationType_BlockRead = 2,
	DTC_DCSOperationType_BlockWrite = 3,
	DTC_DCSOperationType_DoubleRead = 4,
	DTC_DCSOperationType_DoubleWrite = 5,
	DTC_DCSOperationType_Unknown = 0xF
};

/// <summary>
/// Convert a DTC_DCSOperationType enumeration value to its string or JSON representation
/// </summary>
struct DTC_DCSOperationTypeConverter
{
	DTC_DCSOperationType type_;  ///< DTC)DCSOperationType to convert

	/// <summary>
	/// DTC_DCSOperationTypeConverter Constructor
	/// </summary>
	/// <param name="type">DTC_DCSOperationType to convert</param>
	explicit DTC_DCSOperationTypeConverter(DTC_DCSOperationType type)
		: type_(type) {}

	/// <summary>
	/// Convert the type to its string representation
	/// </summary>
	/// <returns></returns>
	std::string toString() const
	{
		switch (type_)
		{
			case DTC_DCSOperationType_Read:
				return "Read";
			case DTC_DCSOperationType_Write:
				return "Write";
			case DTC_DCSOperationType_BlockRead:
				return "BlockRead";
			case DTC_DCSOperationType_BlockWrite:
				return "BlockWrite";
			case DTC_DCSOperationType_DoubleRead:
				return "DoubleRead";
			case DTC_DCSOperationType_DoubleWrite:
				return "DoubleWrite";
			case DTC_DCSOperationType_Unknown:
			default:
				return "Unknown";
		}
	}

	/// <summary>
	/// Write a DTC_DCSOperationTypeConverter in JSON format to the given stream
	/// </summary>
	/// <param name="stream">Stream to write</param>
	/// <param name="type">DTC_DCSOperationTypeConverter to serialize</param>
	/// <returns>Stream reference for continued streaming</returns>
	friend std::ostream& operator<<(std::ostream& stream, const DTC_DCSOperationTypeConverter& type)
	{
		switch (type.type_)
		{
			case DTC_DCSOperationType_Read:
				stream << "\"Read\"";
				break;
			case DTC_DCSOperationType_Write:
				stream << "\"Write\"";
				break;
			case DTC_DCSOperationType_BlockRead:
				stream << "\"BlockRead\"";
				break;
			case DTC_DCSOperationType_BlockWrite:
				stream << "\"BlockWrite\"";
				break;
			case DTC_DCSOperationType_DoubleRead:
				stream << "\"DoubleRead\"";
				break;
			case DTC_DCSOperationType_DoubleWrite:
				stream << "\"DoubleWrite\"";
				break;
			case DTC_DCSOperationType_Unknown:
			default:
				stream << "\"Unknown\"";
				break;
		}
		return stream;
	}
};

/// <summary>
/// The DTC_DataPacket class represents the 16 bytes of raw data for all DTC packets.
/// The class works in two modes: "overlay" mode, where the data is in a fixed location in memory and modification is
/// restricted, and "owner" mode, where the DataPacket is a concrete instance.
/// </summary>
class DTC_DataPacket
{
public:
	/// <summary>
	/// Construct a DTC_DataPacket in owner mode
	/// </summary>
	DTC_DataPacket();

	/// <summary>
	/// Construct a DTC_DataPacket using a pointer to data. Flag will be set that the packet is read-only.
	/// </summary>
	/// <param name="data">Pointer to data</param>
	explicit DTC_DataPacket(const void* data)
		: dataPtr_(static_cast<const uint8_t*>(data)), dataSize_(16), memPacket_(true) {}

	/// <summary>
	/// Creates a copy of the DTC_DataPacket. Mode is preserved, if the existing DataPacket was in "owner" mode, a deep
	/// copy is made, otherwise the reference to the read-only memory will be copied.
	/// </summary>
	/// <param name="in">Input DTC_DataPacket</param>
	DTC_DataPacket(const DTC_DataPacket& in);
	/// <summary>
	/// Default move constructor
	/// </summary>
	/// <param name="in">DTC_DataPacket rvalue</param>
	DTC_DataPacket(DTC_DataPacket&& in) = default;

	virtual ~DTC_DataPacket();

	/// <summary>
	/// Default copy-assignment operator
	/// </summary>
	/// <param name="in">DTC_DataPacket lvalue</param>
	/// <returns>DTC_DataPacket reference</returns>
	DTC_DataPacket& operator=(const DTC_DataPacket& in) = default;
	/// <summary>
	/// Default move-assignment operator
	/// </summary>
	/// <param name="in">DTC_DataPacket rvalue</param>
	/// <returns>DTC_DataPacket reference</returns>
	DTC_DataPacket& operator=(DTC_DataPacket&& in) = default;

	/// <summary>
	/// Set the given word of the DataPacket.
	/// No-op if the DataPacket is in overlay mode
	/// </summary>
	/// <param name="index">Index of the word to change</param>
	/// <param name="data">Value of the word</param>
	void SetWord(uint16_t index, uint8_t data);
	/// <summary>
	/// Gets the current value of the given word
	/// </summary>
	/// <param name="index">Index of the word</param>
	/// <returns>Value of the word</returns>
	uint8_t GetWord(uint16_t index) const;
	/// <summary>
	/// Creates a JSON representation of the DTC_DataPacket
	/// </summary>
	/// <returns>JSON-formatted string representation of the DTC_DataPacket</returns>
	std::string toJSON() const;
	/// <summary>
	/// Create a "packet format" representation of the DTC_DataPacket. See "DTC Hardware User's Guide" for "packet format"
	/// representation.
	/// </summary>
	/// <returns>"packet format" string representation of the DTC_DataPacket</returns>
	std::string toPacketFormat() const;
	/// <summary>
	/// Resize a DTC_DataPacket in "owner" mode. New size must be larger than current.
	/// </summary>
	/// <param name="dmaSize">Size in bytes of the new packet</param>
	/// <returns>If the resize operation was successful</returns>
	bool Resize(const uint16_t dmaSize);

	/// <summary>
	/// Get the current size, in bytes, of the DTC_DataPacket (default: 16)
	/// </summary>
	/// <returns></returns>
	uint16_t GetSize() const { return dataSize_; }

	/// <summary>
	/// Determine whether the DataPacket is in owner mode or overlay mode
	/// </summary>
	/// <returns>True if the DataPacket is in overlay mode</returns>
	bool IsMemoryPacket() const { return memPacket_; }

	/// <summary>
	/// Add a DTC_DataPacket's data to this DataPacket. DataPacket must be large enough to accomodate data and in "owner"
	/// mode.
	/// </summary>
	/// <param name="other">Packet to add</param>
	/// <param name="offset">Where to copy packet in this DataPacket</param>
	/// <returns>True if successful</returns>
	bool CramIn(DTC_DataPacket& other, int offset)
	{
		if (other.dataSize_ + offset <= dataSize_)
		{
			memcpy(const_cast<uint8_t*>(dataPtr_) + offset, other.dataPtr_, other.dataSize_);
			return true;
		}
		return false;
	}

	/// <summary>
	/// Gets the pointer to the data
	/// </summary>
	/// <returns>Pointer to DTC_DataPacket data. Use GetSize() to determine the valid range of this pointer</returns>
	const uint8_t* GetData() const { return dataPtr_; }

	/// <summary>
	/// Comparison operator. Returns this.Equals(other)
	/// </summary>
	/// <param name="other">DataPacket to compare</param>
	/// <returns>this.Equals(other)</returns>
	bool operator==(const DTC_DataPacket& other) const { return Equals(other); }

	/// <summary>
	/// Comparison operator. Returns !this.Equals(other)
	/// </summary>
	/// <param name="other">DataPacket to compare</param>
	/// <returns>!this.Equals(other)</returns>
	bool operator!=(const DTC_DataPacket& other) const { return !Equals(other); }

	/// <summary>
	/// Compare the contents of two DataPackets. Ignores the first two bytes as they may differ (reserved words in most
	/// DMA packets).
	/// </summary>
	/// <param name="other">Data packet to compare</param>
	/// <returns>Whether the two DataPackets contents are equal</returns>
	bool Equals(const DTC_DataPacket& other) const;

	/// <summary>
	/// Serialize a DTC_DataPacket to the given ostream
	/// </summary>
	/// <param name="s">Stream to write to</param>
	/// <param name="p">DataPacket to stream, in binary</param>
	/// <returns>Stream reference for continued streaming</returns>
	friend std::ostream& operator<<(std::ostream& s, DTC_DataPacket& p)
	{
		return s.write(reinterpret_cast<const char*>(p.dataPtr_), p.dataSize_);
	}

private:
	const uint8_t* dataPtr_;
	uint16_t dataSize_;
	bool memPacket_;
	std::vector<uint8_t> vals_;
};

/// <summary>
/// Header information common to all DTC Packets (except Data Packets)
/// </summary>
class DTC_DMAPacket
{
protected:
	uint16_t byteCount_;         ///< Byte count of current block
	bool valid_;                 ///< Whether the DTC believes the packet to be valid
	uint8_t subsystemID_;        ///< Subsystem ID (Data Header packet only)
	DTC_Link_ID linkID_;         ///< Link identifier of packet
	DTC_PacketType packetType_;  ///< Packet type
	uint8_t hopCount_;           ///< Hop count
public:
	/// <summary>
	/// DTC_DMAPacket default constructor. Fills in header fields with default (invalid) values.
	/// </summary>
	DTC_DMAPacket()
		: byteCount_(0), valid_(false), subsystemID_(0), linkID_(DTC_Link_Unused), packetType_(DTC_PacketType_Invalid), hopCount_(0) {}

	/// <summary>
	/// Create a DTC_DMAPacket with the given parameters
	/// </summary>
	/// <param name="type">Packet Type</param>
	/// <param name="link">Link ID</param>
	/// <param name="byteCount">Block byte count. Default is one packet, 16 bytes</param>
	/// <param name="valid">Valid flag for packet, default true</param>
	/// <param name="subsystemID">Subsystem ID for packet</param>
	/// <param name="hopCount">Hop count for packet, default 0</param>
	DTC_DMAPacket(DTC_PacketType type, DTC_Link_ID link, uint16_t byteCount = 16, bool valid = true, uint8_t subsystemID = 0, uint8_t hopCount = 0);

	/// <summary>
	/// Construct a DTC_DMAPacket using the data in the given DataPacket
	/// </summary>
	/// <param name="in">DTC_DataPacket to interpret</param>
	explicit DTC_DMAPacket(const DTC_DataPacket in);
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="in">DTC_DMAPacket to copy</param>
	DTC_DMAPacket(const DTC_DMAPacket& in) = default;
	/// <summary>
	/// Default move Constructor
	/// </summary>
	/// <param name="in">DTC_DMAPacket rvalue</param>
	DTC_DMAPacket(DTC_DMAPacket&& in) = default;

	virtual ~DTC_DMAPacket() = default;

	/// <summary>
	/// Default Copy Assignment Operator
	/// </summary>
	/// <param name="in">DTC_DMAPacket to copy</param>
	/// <returns>DTC_DMAPacket reference</returns>
	DTC_DMAPacket& operator=(const DTC_DMAPacket& in) = default;
	/// <summary>
	/// Default Move Assignment Operator
	/// </summary>
	/// <param name="in">DTC_DMAPacket rvalue</param>
	/// <returns>DTC_DMAPacket reference</returns>
	DTC_DMAPacket& operator=(DTC_DMAPacket&& in) = default;

	/// <summary>
	/// Convert a DTC_DMAPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DMA Header bytes set</returns>
	virtual DTC_DataPacket ConvertToDataPacket() const;

	/// <summary>
	/// Packet Type accessor
	/// </summary>
	/// <returns>Packet Type of DMA Packet</returns>
	DTC_PacketType GetPacketType() const { return packetType_; }

	/// <summary>
	/// Gets the DMA Header in JSON
	/// </summary>
	/// <returns>JSON-formatted string representation of DMA Header information</returns>
	std::string headerJSON() const;
	/// <summary>
	/// Gets the DMA header in "packet format" (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DMA header information</returns>
	std::string headerPacketFormat() const;

	/// <summary>
	/// Gets the block byte count
	/// </summary>
	/// <returns>Block byte count of DMA Header</returns>
	uint16_t GetByteCount() const { return byteCount_; }

	/// <summary>
	/// Gets the Link ID of the packet
	/// </summary>
	/// <returns>The Link ID of the packet</returns>
	DTC_Link_ID GetLinkID() const { return linkID_; }

	/// <summary>
	/// Gets the Hop Count of the packet
	///
	/// </summary>
	/// <returns>The Hop count of the packet</returns>
	uint8_t GetHopCount() const
	{
		return hopCount_;
	}
	/// <summary>
	/// Gets the Subsystem ID of the packet
	///
	/// This method should only be used if the packet is a DataHeader packet
	/// </summary>
	/// <returns>Subsystem ID stored in the packet</returns>
	uint8_t GetSubsystemID() const
	{
		return subsystemID_;
	}

	/// <summary>
	/// Converts the DMA Packet to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DMA packet</returns>
	virtual std::string toPacketFormat();
	/// <summary>
	/// Convert the DMA Packet to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DMA packet</returns>
	virtual std::string toJSON();

	/// <summary>
	/// Stream the JSON representation of the DTC_DMAPacket to the given stream
	/// </summary>
	/// <param name="stream">Stream to write JSON data to</param>
	/// <param name="packet">Packet to stream</param>
	/// <returns>Stream reference for continued streaming</returns>
	friend std::ostream& operator<<(std::ostream& stream, DTC_DMAPacket& packet)
	{
		stream << packet.toJSON();
		return stream;
	}
};

/// <summary>
/// Representation of a DCS Request Packet
/// </summary>
class DTC_DCSRequestPacket : public DTC_DMAPacket
{
public:
	/// <summary>
	/// Default Constructor, zeroes out header fields
	/// </summary>
	DTC_DCSRequestPacket();
	/// <summary>
	/// DCSRequestPacket constructor, using given link and roc
	/// </summary>
	/// <param name="link">Link ID for packet</param>
	DTC_DCSRequestPacket(DTC_Link_ID link);
	/// <summary>
	/// Create a DTC_DCSRequestPacket instance with the given fields filled in
	/// </summary>
	/// <param name="link">Link ID for packet</param>
	/// <param name="type">OpCode of packet</param>
	/// <param name="requestAck">Whether to request acknowledement of this operation</param>
	/// <param name="incrementAddress">Whether to increment the address pointer for block reads/writes</param>
	/// <param name="address">Address of ROC register</param>
	/// <param name="data">Data/wordCount for operation</param>
	/// <param name="address2">Address of ROC register</param>
	/// <param name="data2">Data/wordCount for operation</param>
	DTC_DCSRequestPacket(DTC_Link_ID link, DTC_DCSOperationType type, bool requestAck, bool incrementAddress, uint16_t address,
						 uint16_t data = 0x0, uint16_t address2 = 0x0, uint16_t data2 = 0x0);
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="in">DTC_DCSRequestPacket to copy</param>
	DTC_DCSRequestPacket(const DTC_DCSRequestPacket& in) = default;
	/// <summary>
	/// Default Move Constructor
	/// </summary>
	/// <param name="in">DTC_DCSRequestPacket rvalue</param>
	DTC_DCSRequestPacket(DTC_DCSRequestPacket&& in) = default;
	/// <summary>
	/// Construct a DTC_DCSRequestPacket using the data in the given DataPacket
	/// </summary>
	/// <param name="in">DataPacket to parse</param>
	explicit DTC_DCSRequestPacket(const DTC_DataPacket in);

	/// <summary>
	/// Default Copy Assignment Operator
	/// </summary>
	/// <param name="in">DTC_DCSRequestPacket to copy</param>
	/// <returns>DTC_DCSRequestPacket Reference</returns>
	DTC_DCSRequestPacket& operator=(const DTC_DCSRequestPacket& in) = default;
	/// <summary>
	/// Default Move Assignment Operator
	/// </summary>
	/// <param name="in">DTC_DCSRequestPacket rvalue</param>
	/// <returns>DTC_DCSRequestPacket Reference</returns>
	DTC_DCSRequestPacket& operator=(DTC_DCSRequestPacket&& in) = default;

	virtual ~DTC_DCSRequestPacket() = default;

	/// <summary>
	/// Gets the opcode of the DCS Request Packet
	/// </summary>
	/// <returns>Current Opcode of the DCS Request Packet</returns>
	DTC_DCSOperationType GetType() const { return type_; }

	/// <summary>
	/// Read the double operation bit from the DCS Request packet
	/// </summary>
	/// <returns>Whether the Double operation bit is set</returns>
	bool IsDoubleOp() const { return (type_ & 0x4) != 0; }

	/// <summary>
	/// Read the request acknowledgment bit from the DCS Request Packet
	/// </summary>
	/// <returns>Whether the request acknowledgment bit is set</returns>
	bool RequestsAck() const { return requestAck_; }

	/// <summary>
	/// Read the increment address bit from the DCS Request Packet
	/// </summary>
	/// <returns>Whether the increment address bit is set</returns>
	bool IncrementsAddress() const { return incrementAddress_; }

	/// <summary>
	/// Get the request from the DCS Request Packet
	/// </summary>
	/// <param name="secondOp">Whether to read the second request</param>
	/// <returns>Pair of address, data from the given request</returns>
	std::pair<uint16_t, uint16_t> GetRequest(bool secondOp = false)
	{
		if (!secondOp) return std::make_pair(address1_, data1_);
		return std::make_pair(address2_, data2_);
	}

	/// <summary>
	/// Add a second request to the DCS Request Packet
	/// </summary>
	/// <param name="address">Address of the second request</param>
	/// <param name="data">Data for the second request</param>
	void AddRequest(uint16_t address, uint16_t data = 0x0);
	/// <summary>
	/// Set the block write data
	/// </summary>
	/// <param name="data">Vector of 16-bit words to write to the ROC</param>
	void SetBlockWriteData(std::vector<uint16_t> data)
	{
		blockWriteData_ = data;
		UpdatePacketAndWordCounts();
	}
	/// <summary>
	/// Update the packetCount and data_ (Block word count) words in the DCS request packet, based on the type of operation and blockWriteData
	/// </summary>
	void UpdatePacketAndWordCounts();

	/// <summary>
	/// Sets the opcode of the DCS Request Packet
	/// </summary>
	/// <param name="type">Opcode to set</param>
	/// <param name="reqAck">Whether to request acknowledgment of this operation</param>
	/// <param name="incAddress">Whether to increment the address pointer for block reads/writes</param>
	void SetType(DTC_DCSOperationType type, bool reqAck, bool incAddress)
	{
		requestAck_ = reqAck;
		incrementAddress_ = incAddress;
		type_ = type;
	}

	/// <summary>
	/// Convert a DTC_DCSRequestPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DCS Request Packet contents set</returns>
	DTC_DataPacket ConvertToDataPacket() const override;
	/// <summary>
	/// Convert the DCS Request Packet to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DCS Request packet</returns>
	std::string toJSON() override;
	/// <summary>
	/// Converts the DCS Request Packet to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DCS Request packet</returns>
	std::string toPacketFormat() override;

private:
	DTC_DCSOperationType type_;
	bool requestAck_;
	bool incrementAddress_;
	uint16_t packetCount_;
	uint16_t address1_;
	uint16_t data1_;     ///< Also, blockWriteData0_
	uint16_t address2_;  ///< Also, blockWriteData1_
	uint16_t data2_;     ///< Also, blockWriteData2_
	std::vector<uint16_t> blockWriteData_;
};

/// <summary>
/// The DTC Heartbeat Packet (sometimes referred to as a "Readout Request" packet)
/// </summary>
class DTC_HeartbeatPacket : public DTC_DMAPacket
{
public:
	/// <summary>
	/// Construct a DTC_HeartbeatPacket
	/// </summary>
	/// <param name="link">Destination Link</param>
	explicit DTC_HeartbeatPacket(DTC_Link_ID link);
	/// <summary>
	/// Construct a DTC_HeartbeatPacket
	/// </summary>
	/// <param name="link">Destination Link</param>
	/// <param name="event_tag">Timestamp of request</param>
	/// <param name="eventMode">Debug event mode bytes (Default: nullptr) If not null, must be 6 bytes long</param>
	/// <param name="deliveryRingTDC">TDC value from Delivery Ring</param>
	DTC_HeartbeatPacket(DTC_Link_ID link, DTC_EventWindowTag event_tag, DTC_EventMode eventMode = DTC_EventMode(), uint8_t deliveryRingTDC = 0);
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="right">DTC_HeartbeatPacket to copy</param>
	DTC_HeartbeatPacket(const DTC_HeartbeatPacket& right) = default;
	/// <summary>
	/// Default Move Constructor
	/// </summary>
	/// <param name="right">DTC_HeartbeatPacket to move</param>
	DTC_HeartbeatPacket(DTC_HeartbeatPacket&& right) = default;
	/// <summary>
	/// Construct a DTC_HeartbeatPacket from the given DTC_DataPacket
	/// </summary>
	/// <param name="in">DTC_DataPacket to overlay</param>
	explicit DTC_HeartbeatPacket(const DTC_DataPacket in);

	/// <summary>
	/// Default Destructor
	/// </summary>
	virtual ~DTC_HeartbeatPacket() noexcept = default;

	/// <summary>
	/// Get the DTC_EventWindowTag stored in the HeartbeatPacket
	/// </summary>
	/// <returns>Timestamp of Heartbeat</returns>
	DTC_EventWindowTag GetEventWindowTag() const { return event_tag_; }

	/// <summary>
	/// Get the Mode bytes from the Heartbeat packet
	/// </summary>
	/// <returns>5-byte array containing mode bytes</returns>
	virtual DTC_EventMode GetData() { return eventMode_; }

	/// <summary>
	/// Convert a DTC_HeartbeatPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DTC_HeartbeatPacket contents set</returns>
	DTC_DataPacket ConvertToDataPacket() const override;
	/// <summary>
	/// Convert the DTC_HeartbeatPacket to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DTC_HeartbeatPacket</returns>
	std::string toJSON() override;
	/// <summary>
	/// Converts the DTC_HeartbeatPacket to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DTC_HeartbeatPacket</returns>
	std::string toPacketFormat() override;

private:
	DTC_EventWindowTag event_tag_;
	DTC_EventMode eventMode_;
	uint8_t deliveryRingTDC_;
};

/// <summary>
/// The DTC Data Request Packet
/// </summary>
class DTC_DataRequestPacket : public DTC_DMAPacket
{
public:
	/// <summary>
	/// Construct a DTC_DataRequestPacket
	/// </summary>
	/// <param name="link">Destination Link</param>
	/// <param name="debug">Debug Mode flag (Default: true)</param>
	/// <param name="debugPacketCount">Debug Packet Count (Default: 0)</param>
	/// <param name="type">Debug Type (Default: DTC_DebugType_SpecialSequence</param>
	DTC_DataRequestPacket(DTC_Link_ID link, bool debug = true, uint16_t debugPacketCount = 0,
						  DTC_DebugType type = DTC_DebugType_SpecialSequence);
	/// <summary>
	/// Construct a DTC_DataRequestPacket
	/// </summary>
	/// <param name="link">Destination Link</param>
	/// <param name="event_tag">Timestamp to request data for</param>
	/// <param name="debug">Debug Mode flag (Default: true)</param>
	/// <param name="debugPacketCount">Debug Packet Count (Default: 0)</param>
	/// <param name="type">Debug Type (Default: DTC_DebugType_SpecialSequence</param>
	DTC_DataRequestPacket(DTC_Link_ID link, DTC_EventWindowTag event_tag, bool debug = true, uint16_t debugPacketCount = 0,
						  DTC_DebugType type = DTC_DebugType_SpecialSequence);
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="right">DTC_DataRequestPacket to copy</param>
	DTC_DataRequestPacket(const DTC_DataRequestPacket& right) = default;
	/// <summary>
	/// Default Move Constructor
	/// </summary>
	/// <param name="right">DTC_DataRequestPacket to move</param>
	DTC_DataRequestPacket(DTC_DataRequestPacket&& right) = default;
	/// <summary>
	/// Construct a DTC_DataRequestPacket from the given DTC_DataPacket
	/// </summary>
	/// <param name="in">DTC_DataPacket to overlay</param>
	explicit DTC_DataRequestPacket(const DTC_DataPacket in);

	/// <summary>
	/// Get the value of the debug flag
	/// </summary>
	/// <returns>Debug Flag value</returns>
	bool GetDebug() const { return debug_; }

	/// <summary>
	/// Get the Debug type
	/// </summary>
	/// <returns>DTC_DebugType enumeration value</returns>
	DTC_DebugType GetDebugType() const { return type_; }

	/// <summary>
	/// Get the Debug Packet Count
	/// </summary>
	/// <returns>Number of packets requested by Data Request</returns>
	uint16_t GetDebugPacketCount() const { return debugPacketCount_; }

	/// <summary>
	/// Set the Debug Packet Count
	/// </summary>
	/// <param name="count">Number of packets to request</param>
	void SetDebugPacketCount(uint16_t count);

	/// <summary>
	/// Get the timestamp of the request
	/// </summary>
	/// <returns>DTC_EventWindowTag of reqeust</returns>
	DTC_EventWindowTag GetEventWindowTag() const { return event_tag_; }

	/// <summary>
	/// Convert a DTC_DataRequestPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DTC_DataRequestPacket contents set</returns>
	DTC_DataPacket ConvertToDataPacket() const override;
	/// <summary>
	/// Convert the DTC_DataRequestPacket to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DTC_DataRequestPacket</returns>
	std::string toJSON() override;
	/// <summary>
	/// Converts the DTC_DataRequestPacket to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DTC_DataRequestPacket</returns>
	std::string toPacketFormat() override;

private:
	DTC_EventWindowTag event_tag_;
	bool debug_;
	uint16_t debugPacketCount_;
	DTC_DebugType type_;
};

/// <summary>
/// The DCS Reply Packet
/// </summary>
class DTC_DCSReplyPacket : public DTC_DMAPacket
{
public:
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="right">DTC_DCSReplyPacket to copy</param>
	DTC_DCSReplyPacket(const DTC_DCSReplyPacket& right) = default;
	/// <summary>
	/// Default Move Constructor
	/// </summary>
	/// <param name="right">DTC_DCSReplyPacket to move</param>
	DTC_DCSReplyPacket(DTC_DCSReplyPacket&& right) = default;
	/// <summary>
	/// Construct a DTC_DCSReplyPacket from the given DTC_DataPacket
	/// </summary>
	/// <param name="in">DTC_DataPacket to overlay</param>
	explicit DTC_DCSReplyPacket(const DTC_DataPacket in);

	/// <summary>
	/// Get the DCS Operation Type
	/// </summary>
	/// <returns>DTC_DCSOperationType enumeration value</returns>
	DTC_DCSOperationType GetType() const { return type_; }

	/// <summary>
	/// Read the double operation bit from the DCS reply packet
	/// </summary>
	/// <returns>Whether the double operation bit is set</returns>
	bool IsDoubleOperation() const { return doubleOp_; }

	/// <summary>
	/// Get the "request acknowledgment" bit from the DCS Reply packet
	/// </summary>
	/// <returns></returns>
	bool IsAckRequested() const { return requestAck_; }

	/// <summary>
	/// Check if the DCS Receive FIFO is empty
	/// </summary>
	/// <returns>Value of DCS Receive FIFO Empty flag</returns>
	bool DCSReceiveFIFOEmpty() const { return dcsReceiveFIFOEmpty_; }

	/// <summary>
	/// Get the corrupt flag from the status word
	/// </summary>
	/// <returns>Whether the corrupt bit is set</returns>
	bool ROCIsCorrupt() const { return corruptFlag_; }

	/// <summary>
	/// Get the value packet count field
	/// </summary>
	/// <returns>The number of packets in the block read</returns>
	uint16_t GetBlockPacketCount() const { return packetCount_; }

	/// <summary>
	/// Get the reply from the DCSReplyPacket
	/// </summary>
	/// <param name="secondOp">Whether to read the second operation</param>
	/// <returns>Pair of address, data from the reply packet</returns>
	std::pair<uint16_t, uint16_t> GetReply(bool secondOp = false)
	{
		if (!secondOp) return std::make_pair(address1_, data1_);
		return std::make_pair(address2_, data2_);
	}

	/// <summary>
	/// Get the block read data, if any
	/// </summary>
	/// <returns>Vector of 16-bit words</returns>
	std::vector<uint16_t> GetBlockReadData() const { return blockReadData_; }

	/// <summary>
	/// Convert a DTC_DCSReplyPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DTC_DCSReplyPacket contents set</returns>
	DTC_DataPacket ConvertToDataPacket() const override;
	/// <summary>
	/// Convert the DTC_DCSReplyPacket to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DTC_DCSReplyPacket</returns>
	std::string toJSON() override;
	/// <summary>
	/// Converts the DTC_DCSReplyPacket to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DTC_DCSReplyPacket</returns>
	std::string toPacketFormat() override;

private:
	DTC_DCSOperationType type_;
	bool doubleOp_;
	bool requestAck_;
	bool dcsReceiveFIFOEmpty_;
	bool corruptFlag_;
	uint16_t packetCount_;
	uint16_t address1_;
	uint16_t data1_;
	uint16_t address2_;
	uint16_t data2_;
	std::vector<uint16_t> blockReadData_;
};

/// <summary>
/// The DTC Data Header Packet (A Data Header and its associated Data Packets forms a Data Block)
/// </summary>
class DTC_DataHeaderPacket : public DTC_DMAPacket
{
public:
	/// <summary>
	/// Construct a DTC_DataHeaderPacket
	/// </summary>
	/// <param name="link">Link from which packet came</param>
	/// <param name="packetCount">Number of DTC_DataPackets in Data Block</param>
	/// <param name="status">Status of Data Block</param>
	/// <param name="dtcid">DTC ID from which packet came</param>
	/// <param name="subsystemid">Subsystem ID from which packet came</param>
	/// <param name="packetVersion">Version of data format</param>
	/// <param name="event_tag">Timestamp of Data Packet (Default: DTC_Timetstamp())</param>
	/// <param name="evbMode">EVB Mode byte (Default: 0)</param>
	DTC_DataHeaderPacket(DTC_Link_ID link, uint16_t packetCount, DTC_DataStatus status, uint8_t dtcid, DTC_Subsystem subsystemid,
						 uint8_t packetVersion, DTC_EventWindowTag event_tag = DTC_EventWindowTag(), uint8_t evbMode = 0);
	/// <summary>
	/// Default Copy Constructor
	/// </summary>
	/// <param name="right">DTC_DataHeaderPacket to copy</param>
	DTC_DataHeaderPacket(const DTC_DataHeaderPacket& right) = default;
	/// <summary>
	/// Default Move Constructor
	/// </summary>
	/// <param name="right">DTC_DataHeaderPacket to move</param>
	DTC_DataHeaderPacket(DTC_DataHeaderPacket&& right) = default;
	/// <summary>
	/// Construct a DTC_DataHeaderPacket from the given DTC_DataPacket
	/// </summary>
	/// <param name="in">DTC_DataPacket to overlay</param>
	explicit DTC_DataHeaderPacket(const DTC_DataPacket in);

	/// <summary>
	/// Convert a DTC_DataHeaderPacket to DTC_DataPacket in "owner" mode
	/// </summary>
	/// <returns>DTC_DataPacket with DTC_DataHeaderPacket contents set</returns>
	DTC_DataPacket ConvertToDataPacket() const override;

	/// <summary>
	/// Get the EVB Mode word from the Data Header Packet
	/// </summary>
	/// <returns>EVB Mode of Data Block</returns>
	uint8_t GetEVBMode() const { return evbMode_; }

	/// <summary>
	/// Get the Subsystem ID of the Data Block
	/// </summary>
	/// <returns>DTC_Subsystem enumeration value</returns>
	DTC_Subsystem GetSubsystem() const { return static_cast<DTC_Subsystem>(GetSubsystemID()); }

	/// <summary>
	/// Get the DTC ID of the Data Block
	/// </summary>
	/// <returns>DTC ID of Data Block</returns>
	uint8_t GetID() const { return dtcId_; }

	/// <summary>
	/// Get the number of Data Packets in the Data block
	/// </summary>
	/// <returns>The number of packets in the Data Block</returns>
	uint16_t GetPacketCount() const { return packetCount_; }

	/// <summary>
	/// Get the Data Packet Version identifier from the Data Header
	/// </summary>
	/// <returns>Version number of Data Packets</returns>
	uint8_t GetVersion() const { return dataPacketVersion_; }

	/// <summary>
	/// Get the Timestamp of the Data Block
	/// </summary>
	/// <returns>timestamp of Data Block</returns>
	DTC_EventWindowTag GetEventWindowTag() const { return event_tag_; }

	/// <summary>
	/// Get the Data Status of the Data Block
	/// </summary>
	/// <returns>DTC_DataStatus enumeration value</returns>
	DTC_DataStatus GetStatus() const { return status_; }

	/// <summary>
	/// Convert the DTC_DataHeaderPacket to JSON representation
	/// </summary>
	/// <returns>JSON-formatted string representation of DTC_DataHeaderPacket</returns>
	std::string toJSON() override;
	/// <summary>
	/// Converts the DTC_DataHeaderPacket to "packet format" representation (See DTC_DataPacket::toPacketFormat())
	/// </summary>
	/// <returns>"packet format" string representation of DTC_DataHeaderPacket</returns>
	std::string toPacketFormat() override;

	/// <summary>
	/// Determine if two Data Header packets are equal (Evaluates DataPacket == DataPacket, see DTC_DataPacket::Equals)
	/// </summary>
	/// <param name="other">DataHeaderPacket to compare</param>
	/// <returns>True if Data Header packet contents are equal</returns>
	bool operator==(const DTC_DataHeaderPacket& other) const { return Equals(other); }

	/// <summary>
	/// Determine if two Data Header packets are not equal (Returns !(dhp == dhp))
	/// </summary>
	/// <param name="other">DataHeaderPacket to compare</param>
	/// <returns>True if Data Header packet contents are not equal</returns>
	bool operator!=(const DTC_DataHeaderPacket& other) const { return !Equals(other); }

	/// <summary>
	/// Determine if two Data Header packets are equal (Evaluates DataPacket == DataPacket, see DTC_DataPacket::Equals)
	/// </summary>
	/// <param name="other">DataHeaderPacket to compare</param>
	/// <returns>True if Data Header packet contents are equal</returns>
	bool Equals(const DTC_DataHeaderPacket& other) const;

private:
	uint16_t packetCount_;
	DTC_EventWindowTag event_tag_;
	DTC_DataStatus status_;
	uint8_t dataPacketVersion_;
	uint8_t dtcId_;
	uint8_t evbMode_;
};

/// <summary>
/// A Data Block object (DataHeader packet plus associated Data Packets)
/// Constructed as a pointer to a region of memory
/// </summary>
struct DTC_DataBlock
{
	std::shared_ptr<std::vector<uint8_t>> allocBytes{nullptr};  ///< Used if the block owns its memory
	const void* blockPointer{nullptr};                          ///< Pointer to DataBlock in Memory
	size_t byteSize{0};                                         ///< Size of DataBlock
	mutable std::shared_ptr<DTC_DataHeaderPacket> hdr{nullptr};

	/**
	 * @brief Create a DTC_DataBlock using a pointer to a memory location containing a Data Block
	 * @param ptr Pointer to Data Block
	 *
	 * WARNING: This function assumes that the pointer is pointing to a valid DTC_DataHeaderPacket!
	*/
	DTC_DataBlock(const void* ptr)
		: blockPointer(ptr)
	{
		DTC_DataPacket pkt(ptr);
		DTC_DataHeaderPacket hdr(pkt);
		byteSize = hdr.GetByteCount();
	}

	/// <summary>
	/// Create a DTC_DataBlock pointing to the given location in memory with the given size
	/// </summary>
	/// <param name="ptr">Pointer to DataBlock in memory</param>
	/// <param name="sz">Size of DataBlock</param>
	DTC_DataBlock(const void* ptr, size_t sz)
		: blockPointer(ptr), byteSize(sz) {}

	DTC_DataBlock(size_t sz)
		: allocBytes(new std::vector<uint8_t>(sz)), blockPointer(allocBytes->data()), byteSize(sz)
	{
	}

	inline std::shared_ptr<DTC_DataHeaderPacket> GetHeader() const
	{
		assert(byteSize >= 16);
		if (hdr) return hdr;
		hdr = std::make_shared<DTC_DataHeaderPacket>(DTC_DataPacket(blockPointer));
		return hdr;
	}

	inline const void* GetData() const
	{
		assert(byteSize > 16);
		return static_cast<const void*>(reinterpret_cast<const uint8_t*>(blockPointer) + 16);
	}
};

struct DTC_SubEventHeader
{
	uint64_t inclusive_subevent_byte_count : 25;
	uint64_t reserved1 : 7;
	uint64_t event_tag_low : 32;

	uint64_t event_tag_high : 16;
	uint64_t num_rocs : 8;
	uint64_t event_mode : 40;

	uint64_t dtc_mac : 8;
	uint64_t partition_id : 8;
	uint64_t evb_mode : 8;
	uint64_t source_dtc_id : 8;
	uint64_t reserved2 : 32;

	uint64_t link0_status : 8;
	uint64_t link1_status : 8;
	uint64_t link2_status : 8;
	uint64_t link3_status : 8;
	uint64_t link4_status : 8;
	uint64_t link5_status : 8;
	uint64_t reserved3 : 8;
	uint64_t emtdc : 8;

	DTC_SubEventHeader()
		: inclusive_subevent_byte_count(0)
		, reserved1(0)
		, event_tag_low(0)
		, event_tag_high(0)
		, num_rocs(0)
		, event_mode(0)
		, dtc_mac(0)
		, partition_id(0)
		, evb_mode(0)
		, source_dtc_id(0)
		, reserved2(0)
		, link0_status(0)
		, link1_status(0)
		, link2_status(0)
		, link3_status(0)
		, link4_status(0)
		, link5_status(0)
		, reserved3(0)
		, emtdc(0)
	{}

	std::string toJson() const;
};

class DTC_SubEvent
{
public:
	/// <summary>
	/// Construct a DTC_SubEvent using a pointer to data. Flag will be set that the packet is read-only.
	/// </summary>
	/// <param name="ptr">Pointer to data</param>
	explicit DTC_SubEvent(const uint8_t*& ptr);

	DTC_SubEvent()
		: header_(), data_blocks_() {}

	size_t GetSubEventByteCount() { return header_.inclusive_subevent_byte_count; }

	DTC_EventWindowTag GetEventWindowTag() const;
	void SetEventWindowTag(DTC_EventWindowTag const& tag);
	void SetEventMode(DTC_EventMode const& mode);
	uint8_t GetDTCID() const;

	std::vector<DTC_DataBlock> const& GetDataBlocks() const
	{
		return data_blocks_;
	}
	size_t GetDataBlockCount() const { return data_blocks_.size(); }
	DTC_DataBlock* GetDataBlock(size_t idx)
	{
		if (idx >= data_blocks_.size()) throw std::out_of_range("Index " + std::to_string(idx) + " is out of range (max: " + std::to_string(data_blocks_.size() - 1) + ")");
		return &data_blocks_[idx];
	}
	void AddDataBlock(DTC_DataBlock blk)
	{
		data_blocks_.push_back(blk);
		header_.num_rocs++;
		UpdateHeader();
	}

	DTC_Subsystem GetSubsystem() const { return static_cast<DTC_Subsystem>((header_.source_dtc_id & 0x70) >> 4); }
	void SetSourceDTC(uint8_t id, DTC_Subsystem subsystem = DTC_Subsystem_Other)
	{
		header_.source_dtc_id = (id & 0xf) + ((static_cast<int>(subsystem) & 0x7) << 4);
	}
	DTC_SubEventHeader* GetHeader() { return &header_; }
	void UpdateHeader();

private:
	DTC_SubEventHeader header_;
	std::vector<DTC_DataBlock> data_blocks_;
};

struct DTC_EventHeader
{
	uint64_t inclusive_event_byte_count : 24;
	uint64_t reserved1 : 8;
	uint64_t event_tag_low : 32;

	uint64_t event_tag_high : 16;
	uint64_t num_dtcs : 8;
	uint64_t event_mode : 40;

	uint64_t dtc_mac : 8;
	uint64_t partition_id : 8;
	uint64_t evb_mode : 8;
	uint64_t evb_id : 8;
	uint64_t evb_status : 8;
	uint64_t emtdc : 8;
	uint64_t reserved2 : 16;

	DTC_EventHeader()
		: inclusive_event_byte_count(0)
		, reserved1(0)
		, event_tag_low(0)
		, event_tag_high(0)
		, num_dtcs(0)
		, event_mode(0)
		, dtc_mac(0)
		, partition_id(0)
		, evb_mode(0)
		, evb_id(0)
		, evb_status(0)
		, emtdc(0)
		, reserved2(0)
	{}

	std::string toJson() const;
};

class DTC_Event
{
public:
	/// <summary>
	/// Construct a DTC_Event in "overlay" mode using the given DMA buffer pointer. Flag will be set that the packet
	/// is read-only.
	/// </summary>
	/// <param name="data">Pointer data</param>
	explicit DTC_Event(const void* data);

	explicit DTC_Event(size_t data_size);

	DTC_Event()
		: header_(), sub_events_(), buffer_ptr_(nullptr) {}

	static const int MAX_DMA_SIZE = 0x8000;

	void SetupEvent();
	size_t GetEventByteCount() const { return header_.inclusive_event_byte_count; }
	DTC_EventWindowTag GetEventWindowTag() const;
	void SetEventWindowTag(DTC_EventWindowTag const& tag);
	void SetEventMode(DTC_EventMode const& mode);
	const void* GetRawBufferPointer() const { return buffer_ptr_; }

	std::vector<DTC_SubEvent> const& GetSubEvents() const
	{
		return sub_events_;
	}
	size_t GetSubEventCount() const { return sub_events_.size(); }

	size_t GetSubEventCount(DTC_Subsystem subsys) const
	{
		size_t count = 0;
		for (size_t ii = 0; ii < sub_events_.size(); ++ii)
		{
			if (sub_events_[ii].GetSubsystem() == subsys) ++count;
		}
		return count;
	}

	size_t GetBlockCount(DTC_Subsystem subsys) const
	{
		size_t count = 0;
		for (size_t ii = 0; ii < sub_events_.size(); ++ii)
		{
			if (sub_events_[ii].GetSubsystem() == subsys)
			{
				count += sub_events_[ii].GetDataBlockCount();
			}
		}
		return count;
	}

	DTC_SubEvent* GetSubEvent(size_t idx)
	{
		if (idx >= sub_events_.size()) throw std::out_of_range("Index " + std::to_string(idx) + " is out of range (max: " + std::to_string(sub_events_.size() - 1) + ")");
		return &sub_events_[idx];
	}

	void AddSubEvent(DTC_SubEvent subEvt)
	{
		sub_events_.push_back(subEvt);
		header_.num_dtcs++;
		UpdateHeader();
	}
	DTC_SubEvent* GetSubEventByDTCID(uint8_t dtc, DTC_Subsystem subsys)
	{
		auto dtcid = (dtc & 0xF) + ((static_cast<uint8_t>(subsys) & 0x7) << 4);
		for (size_t ii = 0; ii < sub_events_.size(); ++ii)
		{
			if (sub_events_[ii].GetDTCID() == dtcid) return &sub_events_[ii];
		}
		return nullptr;
	}

	std::vector<DTC_DataBlock> GetSubsystemDataBlocks(DTC_Subsystem subsys) const {
		std::vector<DTC_DataBlock> output;

		for(auto& subevt : sub_events_) {
//FIXME: CRV Subevents don't have the correct subsystem encoded. Skipping the check for now.
//			if(subevt.GetSubsystem() == subsys) {
				auto subevtblocks = subevt.GetDataBlocks();
				output.insert(output.end(), subevtblocks.begin(), subevtblocks.end());
//			}
		}

		return output;
	}

	DTC_EventHeader* GetHeader() { return &header_; }

	void UpdateHeader();
	void WriteEvent(std::ostream& output, bool includeDMAWriteSize = true);

private:
	std::shared_ptr<std::vector<uint8_t>> allocBytes{nullptr};  ///< Used if the block owns its memory
	DTC_EventHeader header_;
	std::vector<DTC_SubEvent> sub_events_;
	const void* buffer_ptr_;
};

}  // namespace DTCLib

#endif  // DTC_PACKETS_H
